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

TFcore::TFcore() {
	M = 1;
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

	free(emgRaw);
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
		xbin[i] = (i / (double)(XRES));
	}

	emgRaw   = (double*)malloc(numch * sizeof(double));
	emgMAV   = (double*)malloc(numch * sizeof(double));
	emgStack = (double*)malloc((int)(WIN_MAV) * numch * sizeof(double));

	xnew = (double*)malloc(numch * sizeof(double));
	idnew = (int*)malloc(numch * sizeof(int));

	tableSize = 500;
	discretizeStep = 100;
	normalTable = (double*)malloc(tableSize * sizeof(double));

	// distribution = max number of patterns
	// realloc vs. malloc?
	pdist_prior = (double*)malloc(XRES * numch * M_MAX);
	pdist_post  = (double*)malloc(XRES * numch * M_MAX);
	pdist_lik   = (double*)malloc(XRES * numch * M_MAX);

	// model parameter
	xmax = 2.00e-4;

	alpha = 1.00e-5;
	beta = 1.00e-10;

	pstar = 1e-17;

	reghold = 512;

	setVal(mu_init, 1.00e-2);
	setVal(sig2_init, 1.00e-3);
	setVal(sig2_reg, 1.00e-2);
	setVal(sig2_update, 1.00e+2);

	setZero(emgRaw);
	setZero(emgMAV);

	constructLookup();

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
	if (isCollect()) {
		bRegister = false;

		reglast += 1;

		pmax = -1;
		mpred = 0;

		for (int n = 0; n < numch; n++) {
			xnew[n] = (1 / xmax) * emgMAV[n];
			xnew[n] = min((double)1, xnew[n]);
			idnew[n] = floor(xnew[n] * (double)(XRES));
			if (idnew[n] > (XRES - 1)) idnew[n] = XRES - 1;
		}

		computeDiffusion();
		computeNormal(xnew, sig2_update, pdist_lik);

		computeProduct<<< XRES, numch * M >>>(pdist_prior, pdist_lik, pdist_post);

		bCollect = false;
		bCompute = true;
	}
	else {
		return;
	}
}

void TFcore::computeDiffusion() {
	// XRES X NUMCH X M
	double* pshift = (double*)malloc(_msize(pdist_prior) * sizeof(double));
	double* dpshift = (double*)malloc(_msize(pdist_prior) * sizeof(double));

	double* dp = (double*)malloc(_msize(pdist_prior) * sizeof(double));
	double* ddp = (double*)malloc(_msize(pdist_prior) * sizeof(double));
	
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numch; n++) {
			pshift[0 + n * XRES + m * (XRES * numch)] = pshift[1 + n * XRES + m * (XRES * numch)];
			memcpy(&pshift[1 + n * XRES + m * (XRES * numch)],
				&pdist_prior[0 + n * XRES + m * (XRES * numch)],
				(XRES - 1) * sizeof(double));
		}
	}
	computeGradient <<< XRES, numch * M >>> (pdist_prior, pshift, dp, (1 / (double)(XRES)));

	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numch; n++) {
			dpshift[0 + n * XRES + m * (XRES * numch)] = dpshift[1 + n * XRES + m * (XRES * numch)];
			memcpy(&dpshift[1 + n * XRES + m * (XRES * numch)],
				&dp[0 + n * XRES + m * (XRES * numch)],
				(XRES - 1) * sizeof(double));
		}
	}
	computeGradient <<< XRES, numch * M >>> (dp, dpshift, ddp, (1 / (double)(XRES)));
	computeWeightedSum <<< XRES, numch * M >>> (pdist_prior, ddp, alpha, beta);
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

void TFcore::computeNormal(double mu, double sig2, double* p) {
	if (_msize(p) == XRES) {
		double Z[XRES];
		
		for (int i = 0; i < XRES; i++) {
			Z[i] = xbin[i] - mu;
			Z[i] = discretizeStep * Z[i] * (1 / sqrt(sig2));
		}
		for (int i = 0; i < XRES; i++) {
			int id = floor(abs(Z[i]));
			if (id < tableSize)
				p[i] = normalTable[id];
			else
				p[i] = EPSILON;
		}

		normalizeProb(p);
	}
	else {
		printf("[error] sizeof @ computeNormal\n");
		return;
	}
}

void TFcore::computeNormal(double* mu, double* sig2, double* p) {
	double* Z = (double*)malloc((int)(XRES) * numch * sizeof(double));
	double* p0 = (double*)malloc((int)(XRES) * numch * sizeof(double));
	if (_msize(p0) == (int)(XRES) * numch) {
		for (int n = 0; n < numch; n++) {
			for (int m = 0; m < XRES; m++) {
				Z[m + n * XRES] = xbin[m] - mu[n];
				Z[m + n * XRES] = discretizeStep * Z[m + n * XRES] *(1 / sqrt(sig2[n]));

				int id = floor(abs(Z[m + n * XRES]));
				if (id < tableSize)
					p0[m + n * XRES] = normalTable[id];
				else
					p0[m + n * XRES] = EPSILON;
			}
		}

		// normalize (XRES X numch)
		normalizeProb(p0);

		// memcpy to (XRES X numch X M)
		for (int m = 0; m < M; m++) {
			memcpy(&p[m], &p0[0], (int)(XRES)*numch * sizeof(double));
		}
	}
	else {
		printf("[error] sizeof @ computeNormal\n");
		return;
	}
}

void TFcore::normalizeProb(double* p) {
	if (_msize(p) == XRES) {
		double psum = 0;
		for (int i = 0; i < XRES; i++) psum += p[i];
		for (int i = 0; i < XRES; i++) p[i] /= psum;
	}
	else if (_msize(p) == (int)(XRES) * numch) {
		for (int n = 0; n < numch; n++) {
			double psum = 0;
			for (int m = 0; m < XRES; m++) psum += p[m + n * XRES];
			for (int m = 0; m < XRES; m++) p[m + n * XRES] /= psum;
		}
	}
}

void TFcore::constructLookup() {
	for (int i = 0; i < _msize(normalTable); i++) {
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
