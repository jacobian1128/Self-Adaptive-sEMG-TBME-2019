#include "TFcore.cuh"

#include <sstream>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//__global__ void computeGradient(double* p, double* pshift, double* dp, double dx) {
//	int i = blockIdx.x * blockDim.x + threadIdx.x;
//	dp[i] = (p[i] - pshift[i]) / (dx);
//}
//
//__global__ void computeWeightedSum(double* p, double* ddp, double alpha, double beta) {
//	int i = blockIdx.x * blockDim.x + threadIdx.x;
//	p[i] += alpha * ddp[i] + beta * ((1 / (double)XRES) - p[i]);
//}
//

//__global__ void computeLikelihood(double* mu, double* sig2, double* lookup, double discretizeStep, double tableSize, double* pdist_lik) {
//	// <<< numChannels, blockDim.x = XRES >>>
//	int i = blockIdx.x * blockDim.x + threadIdx.x;
//	int n = blockIdx.x;
//	int k = threadIdx.x;n
//
//	double Z = 0;
//	double p0 = 0;
//
//	int id = 0;
//
//	Z = (k / (double)XRES)- mu[n];
//	Z = discretizeStep * Z * (1 / sqrt(sig2[n]));
//
//	id = floor(abs(Z));
//	if (id < tableSize) {
//		p0 = lookup[id];
//	}
//	else {
//		p0 = EPSILON;
//	}	
//
//	pdist_lik[i] = p0;
//}

//__global__ void normalizeProb(double* pdist) {
//	// <<< M, numChannels >>>
//	int numChannels = blockDim.x;
//	int m = blockIdx.x;
//	int n = threadIdx.x;
//
//	double psum = 0;
//	for (int k = 0; k < XRES; k++) psum += pdist[k + n * XRES + m * numChannels * XRES];
//	for (int k = 0; k < XRES; k++) pdist[k + n * XRES + m * numChannels * XRES] /= psum;
//}

__global__ void _computeProduct(double* pdist_prior, double* pdist_lik, double* pdist_post) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	pdist_post[i] = pdist_prior[i] * pdist_lik[i];
}

TFcore::TFcore(CTFcoreMFCDlg& dlg) : dlg_(dlg) {
	M = 1;
	reglast = 0;

	bCollect = false;
	bCompute = false;
	bRegister = false;

	bGPU_product = false;
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

	outFile.close();
}

void TFcore::initModel(int n) {
	// set channel
	setNumChannels(n);

	// memory allocation
	mu_init = (double*)calloc(numChannels, sizeof(double));
	sig2_init = (double*)calloc(numChannels, sizeof(double));

	sig2_reg = (double*)calloc(numChannels, sizeof(double));
	sig2_update = (double*)calloc(numChannels, sizeof(double));

	emgMAV = (double*)calloc(numChannels, sizeof(double));

	xnew = (double*)calloc(numChannels, sizeof(double));
	idnew = (int*)calloc(numChannels, sizeof(int));

	xbin = (double*)calloc(XRES, sizeof(double));
	for (int k = 0; k < XRES; k++) xbin[k] = ((k + 1) / (double)(XRES));
	emgStack = (double*)calloc(WIN_MAV * numChannels, sizeof(double));

	// lookup table for normal distribution
	tableSize = 50000;
	discretizeStep = 10000;
	normalTable = (double*)calloc(tableSize, sizeof(double));

	constructLookup();

	// distribution = max number of patterns
	pdist_prior = (double*)calloc(XRES * numChannels * M_MAX, sizeof(double));
	pdist_post = (double*)calloc(XRES * numChannels * M_MAX, sizeof(double));
	pdist_lik = (double*)calloc(XRES * numChannels * M_MAX, sizeof(double));

	dp = (double*)calloc(XRES * numChannels * M_MAX, sizeof(double));
	ddp = (double*)calloc(XRES * numChannels * M_MAX, sizeof(double));
	p_shift = (double*)calloc(XRES * numChannels * M_MAX, sizeof(double));
	dp_shift = (double*)calloc(XRES * numChannels * M_MAX, sizeof(double));

	p_lik = (double*)calloc(M_MAX, sizeof(double));

	// model parameter
	xmax = 1.00e+1;

	alpha = 1.00e-5;
	beta = 1.00e-10;

	p_star = 1e-14;

	reghold = 16;

	setVal(mu_init, 2.00e-2);
	setVal(sig2_init, 1.00e-2);
	setVal(sig2_reg, 1.00e-3);
	setVal(sig2_update, 1.00e+1);

	setZero(emgMAV);

	// initialize first pattern
	for (int n = 0; n < numChannels; n++) {
		for (int k = 0; k < XRES; k++) {
			// assign uniform distribution
			pdist_prior[k + n * XRES] = (1 / (double)XRES);
			pdist_post[k + n * XRES] = (1 / (double)XRES);
			pdist_lik[k + n * XRES] = (1 / (double)XRES);
		}
	}

	outFile.open("output.csv", ios::out);
}

void TFcore::getSample() {

	// a vertical stack of row samples to compute MAV
	// the first row = the latest sample
	for (int n = 0; n < numChannels; n++) {
		memcpy(&emgStack[1 + n * WIN_MAV], &emgStack[0 + n * WIN_MAV], (WIN_MAV - 1) * sizeof(double));
		emgStack[n * WIN_MAV] = dataStack.front();
		dataStack.erase(dataStack.begin());
	}

	// compute MAV = a single row vector
	for (int n = 0; n < numChannels; n++) {
		emgMAV[n] = 0;
		for (int m = 0; m < WIN_MAV; m++) {
			emgMAV[n] += emgStack[m + n * WIN_MAV] / (double)(WIN_MAV);
		}
	}

	bCollect = true;
	bCompute = false;
}

void TFcore::proceedIteration() {
	if (isCollect()) {
		bRegister = false;

		reglast += 1;

		p_max = -1;
		mpred = 0;

		for (int n = 0; n < numChannels; n++) {
			xnew[n] = (1 / xmax) * emgMAV[n];
			xnew[n] = min((double)1, xnew[n]);
			idnew[n] = floor(xnew[n] * (double)(XRES));
			if (idnew[n] > (XRES - 1)) idnew[n] = XRES - 1;
		}

		computeDiffusion();
		normalizePrior();

		computeLikelihood(xnew, sig2_update);
		normalizeLikelihood();

		if (bGPU_product) {
			//for (int m = 0; m < M; m++) {
			//	double *dp1, *dp2, *dp3;
			//	cudaMalloc((void**)&dp1, XRES * numChannels * sizeof(double));
			//	cudaMalloc((void**)&dp2, XRES * numChannels * sizeof(double));
			//	cudaMalloc((void**)&dp3, XRES * numChannels * sizeof(double));
			//	cudaMemcpy(dp1, &pdist_prior[m * XRES * numChannels], XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);
			//	cudaMemcpy(dp2, &pdist_lik[m * XRES * numChannels], XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);
			//	_computeProduct <<< XRES, numChannels >>> (dp1, dp2, dp3);
			//	cudaMemcpy(&pdist_post[m * XRES * numChannels], dp3, XRES * numChannels * sizeof(double), cudaMemcpyDeviceToHost);
			//	cudaFree(dp1);
			//	cudaFree(dp2);
			//	cudaFree(dp3);
			//}

			double* dp1, * dp2, * dp3;
			cudaMalloc((void**)&dp1, M * XRES * numChannels * sizeof(double));
			cudaMalloc((void**)&dp2, M * XRES * numChannels * sizeof(double));
			cudaMalloc((void**)&dp3, M * XRES * numChannels * sizeof(double));
			cudaMemcpy(dp1, pdist_prior, M * XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dp2, pdist_lik, M * XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);
			_computeProduct <<< M, XRES * numChannels >>> (dp1, dp2, dp3);
			cudaMemcpy(pdist_post, dp3, M * XRES * numChannels * sizeof(double), cudaMemcpyDeviceToHost);
			cudaFree(dp1);
			cudaFree(dp2);
			cudaFree(dp3);
		}
		else {
			computeProduct();
		}
		normalizePost();

		// get maximum likelihood probability
		for (int m = 0; m < M; m++) {
			p_lik[m] = 1;
			for (int n = 0; n < numChannels; n++) {
				p_lik[m] *= pdist_post[idnew[n] + n * XRES + m * numChannels * XRES];
			}
			if (p_lik[m] < 1e-200) p_lik[m] = 1e-200;

			// MLE prediction
			if (p_lik[m] > p_max) {
				p_max = p_lik[m];
				mpred = m;
			}
		}

		// registration
		if ((p_max < p_star) && (reglast > reghold)) {
			registerPattern(xnew, sig2_reg);
			p_lik[M] = 1;
			for (int n = 0; n < numChannels; n++) {
				p_lik[M] *= pdist_post[idnew[n] + n * XRES + M * numChannels * XRES];
			}

			mpred = M;
			M += 1;
			reglast = 0;

			bRegister = true;
		}

		memcpy(&pdist_prior[mpred * XRES * numChannels], &pdist_post[mpred * XRES * numChannels], XRES * numChannels * sizeof(double));

		bCollect = false;
		bCompute = true;
	}
	else {
		return;
	}
}

void TFcore::writeResult() {
	outFile << p_max << ',' << mpred << ',';
	for (int n = 0; n < numChannels; n++) outFile << xnew[n] << ',';
	for (int n = 0; n < numChannels; n++) outFile << emgMAV[n] << ',';
	for (int k = 0; k < XRES; k++) outFile << pdist_prior[k] << ',';
	for (int k = 0; k < XRES; k++) outFile << pdist_lik[k] << ',';
	for (int k = 0; k < XRES; k++) outFile << pdist_post[k] << ',';
	outFile << numSamples << ',' << numChannels;
	outFile << endl;
}

void TFcore::computeDiffusion() {
	// memcpy pdist_prior => pshift(dist_prior)
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			p_shift[0 + n * XRES + m * (XRES * numChannels)] = p_shift[1 + n * XRES + m * (XRES * numChannels)];
			memcpy(&p_shift[1 + n * XRES + m * (XRES * numChannels)],
				&pdist_prior[0 + n * XRES + m * (XRES * numChannels)],
				(XRES - 1) * sizeof(double));
		}
	}
	// compute (p, pshift) => dp
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			dp[0 + n * XRES + m * XRES * numChannels] = 0;
			for (int k = 1; k < XRES; k++) {
				dp[k + n * XRES + m * XRES * numChannels]
					= (1 / XRES) * (pdist_prior[k + n * XRES + m * XRES * numChannels] - p_shift[k - 1 + n * XRES + m * XRES * numChannels]);
			}
		}
	}


	// memcpy dp(dist_prior) => dpshift(dist_prior)
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			dp_shift[0 + n * XRES + m * (XRES * numChannels)] = dp_shift[1 + n * XRES + m * (XRES * numChannels)];
			memcpy(&dp_shift[1 + n * XRES + m * (XRES * numChannels)],
				&dp[0 + n * XRES + m * (XRES * numChannels)],
				(XRES - 1) * sizeof(double));
		}
	}

	// compute (dp, dpshift) => ddp
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			ddp[0 + n * XRES + m * XRES * numChannels] = 0;
			for (int k = 1; k < XRES; k++) {
				ddp[k + n * XRES + m * XRES * numChannels]
					= (1 / XRES) * (dp[k + n * XRES + m * XRES * numChannels] - dp_shift[k + n * XRES + m * XRES * numChannels]);
			}
		}
	}

	// weighted summation
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			for (int k = 0; k < XRES; k++) {
				pdist_prior[k + n * XRES + m * XRES * numChannels]
					= pdist_prior[k + n * XRES + m * XRES * numChannels] + alpha * ddp[k + n * XRES + m * XRES * numChannels];
			}
		}
	}

}

void TFcore::computeProduct() {
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			for (int k = 0; k < XRES; k++) {
				pdist_post[k + n * XRES + m * XRES * numChannels]
					= pdist_prior[k + n * XRES + m * XRES * numChannels]
					* pdist_lik[k + n * XRES + m * XRES * numChannels];
			}
		}
	}
}

void TFcore::computeLikelihood(double* mu, double* sig2) {
	double* Z = (double*)calloc(XRES * numChannels, sizeof(double));
	int id = 0;
	for (int k = 0; k < XRES; k++) {
		for (int n = 0; n < numChannels; n++) {
			Z[k + n * XRES] = ((k + 1) / (double)XRES) - mu[n];
			Z[k + n * XRES] *= discretizeStep / sqrt(sig2[n]);

			id = floor(abs(Z[k + n * XRES]));
			if (id < tableSize) {
				pdist_lik[k + n * XRES] = normalTable[id];
			}
			else {
				pdist_lik[k + n * XRES] = EPSILON;
			}
		}
	}

	for (int m = 1; m < M; m++) memcpy(&pdist_lik[m], &pdist_lik[0], XRES * numChannels * sizeof(double));
	free(Z);
}

void TFcore::loadData(string filename) {
	string extension = filename.substr(filename.find_last_of(".") + 1);
	//if (!strcmp(extension.c_str(), "mat")) {
	//	MATFile* pmat = matOpen(filename.c_str(), "r");

	//	if (pmat == NULL) {
	//		return;
	//	}

	//	mxArray* mxdata = matGetVariable(pmat, "emg");

	//	int m, n;
	//	m = mxGetM(mxdata);
	//	n = mxGetN(mxdata);


	//	double* data = mxGetPr(mxdata);

	//	for (int i = 0; i < m * n; i++) {
	//		dataStack.push_back(data[i]);
	//	}

	//	initModel(n);
	//	setNumSamples(m);

	//	mxDestroyArray(mxdata);
	//	matClose(pmat);
	//}
	if (!strcmp(extension.c_str(), "csv")) {
		// comma separated file (samples x channels)
		fstream file;
		file.open(filename.c_str(), ios::in);

		int m = 0;
		int n = 0;

		string data;
		string line;
		stringstream parse;

		while (getline(file, line, '\n'))
		{
			m++;
			parse.clear();
			parse.str(line);
			while (getline(parse, data, ',')) dataStack.push_back(stod(data));
		}

		n = dataStack.size() / m;

		initModel(n);
		setNumSamples(m);

		file.close();
	}
}

void TFcore::loadData(vector<double> data, int ch) {
	initModel(ch);
	setNumSamples(data.size() / ch);
	for (int k = 0; k < data.size(); k++) dataStack.push_back(data[k]);
}

void TFcore::registerPattern(double* mu, double* sig2) {
	double* Z = (double*)calloc(XRES * numChannels, sizeof(double));
	for (int n = 0; n < numChannels; n++) {
		for (int k = 0; k < XRES; k++) {
			Z[k + n * XRES] = xbin[k] - mu[n];
			Z[k + n * XRES] = discretizeStep * Z[k + n * XRES] * (1 / sqrt(sig2[n]));

			int id = floor(abs(Z[k + n * XRES]));
			if (id < tableSize)
				pdist_post[k + n * XRES + M * XRES * numChannels] = normalTable[id];
			else
				pdist_post[k + n * XRES + M * XRES * numChannels] = EPSILON;
		}
	}

	double psum;
	for (int n = 0; n < numChannels; n++) {
		psum = 0;
		for (int k = 0; k < XRES; k++) psum += pdist_post[k + n * XRES + M * XRES * numChannels];
		for (int k = 0; k < XRES; k++) pdist_post[k + n * XRES + M * XRES * numChannels] /= psum;
	}

	free(Z);
}

void TFcore::normalizeLikelihood() {
	double psum;
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			psum = 0;
			for (int k = 0; k < XRES; k++) psum += pdist_lik[k + n * XRES + m * XRES * numChannels];
			for (int k = 0; k < XRES; k++) pdist_lik[k + n * XRES + m * XRES * numChannels] /= psum;
		}
	}
}

void TFcore::normalizePrior() {
	double psum;
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			psum = 0;
			for (int k = 0; k < XRES; k++) psum += pdist_prior[k + n * XRES + m * XRES * numChannels];
			for (int k = 0; k < XRES; k++) pdist_prior[k + n * XRES + m * XRES * numChannels] /= psum;
		}
	}
}

void TFcore::normalizePost() {
	double psum;
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			psum = 0;
			for (int k = 0; k < XRES; k++) psum += pdist_post[k + n * XRES + m * XRES * numChannels];
			for (int k = 0; k < XRES; k++) pdist_post[k + n * XRES + m * XRES * numChannels] /= psum;
		}
	}
}

void TFcore::constructLookup() {
	for (int i = 0; i < tableSize; i++) {
		normalTable[i] = (1 / sqrt(2 * M_PI)) * exp(-0.5 * ((double)i / discretizeStep) * ((double)i / discretizeStep));
	}
}

void TFcore::setVal(double* vector, double val) {
	for (int k = 0; k < numChannels; k++) {
		vector[k] = val;
	}
}

void TFcore::setZero(double* vector) {
	setVal(vector, (double)0);
}

void TFcore::setOnes(double* vector) {
	setVal(vector, (double)1);
}
