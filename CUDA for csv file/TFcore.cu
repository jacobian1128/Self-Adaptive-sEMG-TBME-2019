#include <mat.h> // for MATLAB .mat file

#include "TFcore.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void computeGradient(double* p, double* pshift, double* dp, double dx) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	dp[i] = (p[i] - pshift[i]) / (dx);
}

__global__ void computeWeightedSum(double* p, double* ddp, double alpha, double beta) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	p[i] += alpha * ddp[i] + beta * ((1 / (double)XRES) - p[i]);
}

__global__ void computeProduct(double* pdist_prior, double* pdist_lik, double* pdist_post) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	pdist_post[i] = pdist_prior[i] * pdist_lik[i];
}

__global__ void computeLikelihood(double* mu, double* sig2, double* lookup, double discretizeStep, double tableSize, double* pdist_lik) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int n = blockIdx.x;
	int k = threadIdx.x;

	double Z = 0;
	double p0 = 0;

	int id = 0;

	Z = (k / (double)XRES)- mu[n];
	Z = discretizeStep * Z * (1 / sqrt(sig2[n]));

	id = floor(abs(Z));
	if (id < tableSize) {
		p0 = lookup[id];
	}
	else {
		p0 = EPSILON;
	}	
	
	pdist_lik[i] = p0;
}

TFcore::TFcore() {
	M = 2000;
	reglast = 0;

	bCollect = false;
	bCompute = false;
	bRegister = false;
}

TFcore::~TFcore() {
	free(mu_init);
	free(sig2_init);
	free(sig2_reg);
	free(sig2_update);

	free(xbin);
	free(xnew);
	free(idnew);

	free(emgMAV);
	free(emgStack);

	free(pdist_prior);
	free(pdist_post);
	free(pdist_lik);
}

void TFcore::initModel(int n) {
	// set channel
	numch = n;

	// memory allocation
	mu_init = (double*)malloc(numch * sizeof(double));
	sig2_init = (double*)malloc(numch * sizeof(double));

	sig2_reg = (double*)malloc(numch * sizeof(double));
	sig2_update = (double*)malloc(numch * sizeof(double));

	xbin = (double*)malloc(XRES * sizeof(double));
	for (int i = 0; i < XRES; i++) {
		xbin[i] = ((i + 1) / (double)(XRES));
	}

	emgMAV   = (double*)malloc(numch * sizeof(double));
	//emgStack = (double*)malloc((int)WIN_MAV * numch * sizeof(double));
	emgStack = (double*)calloc(WIN_MAV * numch, sizeof(double));

	xnew = (double*)malloc(numch * sizeof(double));
	idnew = (int*)malloc(numch * sizeof(int));

	tableSize = 5000;
	discretizeStep = 1000;
	normalTable = (double*)malloc(tableSize * sizeof(double));

	// distribution = max number of patterns
	// calloc vs. malloc?
	pdist_prior = (double*)malloc(XRES * numch * M_MAX * sizeof(double));
	pdist_post  = (double*)malloc(XRES * numch * M_MAX * sizeof(double));
	pdist_lik   = (double*)malloc(XRES * numch * M_MAX * sizeof(double));

	p_lik = (double*)malloc(M_MAX * sizeof(double));

	// model parameter
	xmax = 2.00e-4;

	alpha = 1.00e-5;
	beta = 1.00e-10;

	p_star = 1e-17;

	reghold = 512;

	setVal(mu_init, 1.00e-2);
	setVal(sig2_init, 1.00e-3);
	setVal(sig2_reg, 1.00e-2);
	setVal(sig2_update, 1.00e+2);

	setZero(emgMAV);

	constructLookup();
	printf("initModel done");
}

void TFcore::getSample() {

	// a vertical stack of row samples to compute MAV
	// the first row = the latest sample
	for (int n = 0; n < numch; n++) {
		memcpy(&emgStack[1 + n * WIN_MAV], &emgStack[0 + n * WIN_MAV], (WIN_MAV - 1) * sizeof(double));
		emgStack[n * WIN_MAV] = dataStack.front();
		dataStack.pop_front();
	}
	
	// compute MAV = a single row vector
	for (int n = 0; n < numch; n++) {
		emgMAV[n] = 0;
		for (int m = 0; m < WIN_MAV; m++) {
			emgMAV[n] += emgStack[m + n * WIN_MAV] / (double)(WIN_MAV);
		}
	}
	
	bCollect = true;
	bCompute = false;
}

void TFcore::proceed() {
	if (1) {
	//if (isCollect()) {
		bRegister = false;

		reglast += 1;

		p_max = -1;
		mpred = 0;

		for (int n = 0; n < numch; n++) {
			xnew[n] = (1 / xmax) * emgMAV[n];
			xnew[n] = min((double)1, xnew[n]);
			idnew[n] = floor(xnew[n] * (double)(XRES));
			if (idnew[n] > (XRES - 1)) idnew[n] = XRES - 1;
		}

		printf("diffusion >> ");
		computeDiffusion();
		printf("likelihood >> ");
		//computeNormal(xnew, sig2_update, pdist_lik);
		computeLikelihood <<< numch, XRES >>> (xnew, sig2_update, normalTable, discretizeStep, tableSize, pdist_lik);
		for (int m = 1; m < M; m++) memcpy(&pdist_lik[m], &pdist_lik[0], (int)XRES * numch * sizeof(double));

		printf("posterior >> ");
		computeProduct <<< XRES, numch * M >>> (pdist_prior, pdist_lik, pdist_post);
		normalizeProb(pdist_post);

		printf("prediction >> ");
		for (int m = 0; m < M; m++) {
			p_lik[m] = 1;
			for (int n = 0; n < numch; n++) {
				p_lik[m] *= pdist_post[idnew[n] + n * numch + m * numch * XRES];
			}
			if (p_lik[m] < 1e-99) p_lik[m] = 1e-99;

			// MLE prediction
			if (p_lik[m] > p_max) {
				p_max = p_lik[m];
				mpred = m;
			}
		}


		// registration
		if ((p_max < p_star) && (reglast > reghold)) {
			registerPattern(xnew, sig2_reg, &pdist_post[M * numch * XRES]);
			pdist_post[M] = pdist_prior[M];
			p_lik[M] = 1;
			for (int n = 0; n < numch; n++) {
				p_lik[M] *= pdist_post[idnew[n] + n * numch + M * numch * XRES];
			}

			mpred = M;
			M += 1;
			reglast = 0;

			bRegister = true;
		}

		bCollect = false;
		bCompute = true;
	}
	else {
		//return;
	}
}

void TFcore::computeDiffusion() {
	// XRES X NUMCH X M
	double* pshift = (double*)malloc(XRES * numch * M * sizeof(double));
	double* dpshift = (double*)malloc(XRES * numch * M * sizeof(double));

	double* dp = (double*)malloc(XRES * numch * M * sizeof(double));
	double* ddp = (double*)malloc(XRES * numch * M * sizeof(double));
	
	// memcpy pdist_prior => pshift(dist_prior)
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numch; n++) {
			pshift[0 + n * XRES + m * (XRES * numch)] = pshift[1 + n * XRES + m * (XRES * numch)];
			memcpy(&pshift[1 + n * XRES + m * (XRES * numch)],
				&pdist_prior[0 + n * XRES + m * (XRES * numch)],
				(XRES - 1) * sizeof(double));
		}
	}
	// compute (p, pshift) => dp
	computeGradient <<< 1, XRES * numch * M >>> (pdist_prior, pshift, dp, (1 / (double)(XRES)));

	// memcpy dp(dist_prior) => dpshift(dist_prior)
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numch; n++) {
			dpshift[0 + n * XRES + m * (XRES * numch)] = dpshift[1 + n * XRES + m * (XRES * numch)];
			memcpy(&dpshift[1 + n * XRES + m * (XRES * numch)],
				&dp[0 + n * XRES + m * (XRES * numch)],
				(XRES - 1) * sizeof(double));
		}
	}

	// compute (dp, dpshift) => ddp
	computeGradient <<< 1, XRES* numch* M >>> (dp, dpshift, ddp, (1 / (double)(XRES)));
	computeWeightedSum <<< 1, XRES * numch * M >>> (pdist_prior, ddp, alpha, beta);

	normalizeProb(pdist_prior);
}

void TFcore::loadData(string filename) {
	string extension = filename.substr(filename.find_last_of(".") + 1);
	printf("load %s with extension %s\n", filename.c_str(), extension.c_str());
	if (!strcmp(extension.c_str(), "mat")) {
		printf("load .mat file... ");
		MATFile* pmat = matOpen(filename.c_str(), "r");

		if (pmat == NULL) {
			printf("Error opening file %s\n", filename.c_str());
			return;
		}

		mxArray* mxdata = matGetVariable(pmat, "emg");

		int m, n;
		m = mxGetM(mxdata);
		n = mxGetN(mxdata);

		printf("array size = %d x %d... ", m, n);

		double* data = mxGetPr(mxdata);

		for (int i = 0; i < m * n; i++) {
			dataStack.push_back(data[i]);
		}
		printf("done!\n");

		printf("Memory allocation... ");
		initModel(n);
		setSample(m);
		printf("done!\n");

		mxDestroyArray(mxdata);
		matClose(pmat);
	}
	else if (!strcmp(extension.c_str(), "csv")) {
		printf("load .csv file... ");

		fstream file;
		file.open(filename.c_str(), ios::in);

		int m = 0;
		int n = 0;

		string data;
		while (getline(file, data, '\n')) m++;

		// move to beginning
		file.clear();
		file.seekg(0, ios::beg);

		while (getline(file, data, ','))
		{
			dataStack.push_back(stod(data));
		}

		n = dataStack.size() / m;
		printf("array size = %d x %d... ", m, n);
		printf("done!\n");

		printf("Memory allocation... ");
		initModel(n);
		setSample(m);
		printf("done!\n");

		file.close();
	}
}

void TFcore::registerPattern(double* mu, double* sig2, double* p) {
	double* Z = (double*)malloc((int)XRES * numch * sizeof(double));
	for (int n = 0; n < numch; n++) {
		for (int m = 0; m < XRES; m++) {
			Z[m + n * XRES] = xbin[m] - mu[n];
			Z[m + n * XRES] = discretizeStep * Z[m + n * XRES] * (1 / sqrt(sig2[n]));

			int id = floor(abs(Z[m + n * XRES]));
			if (id < tableSize)
				p[m + n * XRES] = normalTable[id];
			else
				p[m + n * XRES] = EPSILON;
		}
	}

	// normalize (XRES X numch)
	normalizeProb(p);
}

void TFcore::computeNormal(double* mu, double* sig2, double* p) {
	double* Z = (double*)malloc(XRES * numch * sizeof(double));
	double* p0 = (double*)malloc(XRES * numch * sizeof(double));
	
	int id = 0;
	
	//if (_msize(p0) == (int)(XRES) * numch) {
	if (1) {
		for (int n = 0; n < numch; n++) {
			for (int k = 0; k < XRES; k++) {
				Z[k + n * XRES] = xbin[k] - mu[n];
				Z[k + n * XRES] = discretizeStep * Z[k + n * XRES] *(1 / sqrt(sig2[n]));

				id = floor(abs(Z[k + n * XRES]));
				if (id < tableSize) {
					p0[k + n * XRES] = normalTable[id];
				}
				else {
					p0[k + n * XRES] = EPSILON;
				}
			}
		}

		// normalize (XRES X numch)
		//normalizeProb(p0);

		printf("memcpy >> ");
		// memcpy to (XRES X numch X M)
		for (int m = 0; m < M; m++) {
			memcpy(&p[m], &p0[0], (int)XRES * numch * sizeof(double));
		}
	}
	else {
		printf("[error] sizeof @ computeNormal\n");
		//return;
	}

	free(Z);
	free(p0);
}

void TFcore::normalizeProb(double* p) {
	if (_msize(p) == XRES) {
		double psum = 0;
		for (int i = 0; i < XRES; i++) psum += p[i];
		for (int i = 0; i < XRES; i++) p[i] /= psum;
	}
	else if (_msize(p) == (int)XRES * numch) {
		for (int n = 0; n < numch; n++) {
			double psum = 0;
			for (int m = 0; m < XRES; m++) psum += p[m + n * XRES];
			for (int m = 0; m < XRES; m++) p[m + n * XRES] /= psum;
		}
	}
	else if (_msize(p) == (int)XRES * numch * M_MAX) {
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < numch; n++) {
				double psum = 0;
				for (int k = 0; k < XRES; k++) {
					psum += p[k + n * XRES + numch * XRES * m];
				}
				for (int k = 0; k < XRES; k++) {
					p[k + n * XRES + numch * XRES * m] /= psum;
				}
			}
		}
	}
}

void TFcore::constructLookup() {
	for (int i = 0; i < tableSize; i++) {
		normalTable[i] = (1 / sqrt(2 * M_PI)) * exp(-0.5 * ((double)i / discretizeStep) * ((double)i / discretizeStep));
	}
}

void TFcore::setSample(int k) {
	numsample = k;
}

void TFcore::setVal(double* vector, double val) {
	for (int k = 0; k < numch; k++) {
		vector[k] = val;
	}
}

void TFcore::setZero(double* vector) {
	setVal(vector, (double)0);
}

void TFcore::setOnes(double* vector) {
	setVal(vector, (double)1);
}
