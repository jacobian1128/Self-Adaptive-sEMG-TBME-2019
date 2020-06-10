#include "TFcore.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void _computeGradient(double* p, double* pshift, double* dp, int XRES) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	dp[i] = (p[i] - pshift[i]) * (double)XRES;
}

__global__ void _computeWeightedSum(double* p, double* ddp, double alpha) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	//p[i] += alpha * ddp[i] + beta * ((1 / (double)XRES) - p[i]);
	p[i] = p[i] + alpha * ddp[i];
}

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

TFcore::TFcore() {
}

TFcore::TFcore(string filename) {
	outputFile.open(filename, ios::out);
}

TFcore::TFcore(fs::path filename) {
	outputFile.open(filename, ios::out);
}

//TFcore::TFcore(CTFcoreMFCDlg& dlg) : dlg_(dlg) {
//	TFcore();
//}

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

	inputFile.close();
	outputFile.close();
}

void TFcore::initModel(int n) {
	M = 1;
	reglast = 0;
	samples = 0;

	bCollect = false;
	bCompute = false;
	bRegister = false;

	bGPU_product = false;
	bGPU_diffusion = false;

	bFuncRegistration = true;
	bFuncAdaptation = true;

	XRES = 128;
	WIN_MAV = 1024;
	M_MAX = 2048;

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
	tableSize = 10000;
	discretizeStep = 1000;
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
	for (int m = 0; m < M_MAX; m++) {
		p_lik[m] = NAN;
		for (int n = 0; n < numChannels; n++) {
			for (int k = 0; k < XRES; k++) {
				pdist_prior[k + n * XRES + m * XRES * numChannels] = NAN;
				pdist_post[k + n * XRES + m * XRES * numChannels] = NAN;
				pdist_lik[k + n * XRES + m * XRES * numChannels] = NAN;

				dp[k + n * XRES + m * XRES * numChannels] = 0;
				ddp[k + n * XRES + m * XRES * numChannels] = 0;
				p_shift[k + n * XRES + m * XRES * numChannels] = (1 / (double)XRES);
				dp_shift[k + n * XRES + m * XRES * numChannels] = 0;
			}
		}
	}

	// model parameter
	xmax = 5.00e-5;

	alpha = 1.00e-10;
	beta = 1.00e-50;

	p_star = -20;

	reghold = 512;

	setVal(mu_init, 3.00e-2);
	setVal(sig2_init, 1.00e-1);
	setVal(sig2_reg, 1.00e-2);
	setVal(sig2_update, 1.00e+2);

	setZero(emgMAV);

	// initialize first pattern
	for (int n = 0; n < numChannels; n++) {
		for (int k = 0; k < XRES; k++) {
			// assign uniform distribution
			pdist_prior[k + n * XRES] = (1 / (double)XRES);
			pdist_post[k + n * XRES] = (1 / (double)XRES);
			pdist_lik[k + n * XRES] = (1 / (double)XRES);

			dp[k + n * XRES] = 0;
			ddp[k + n * XRES] = 0;
			p_shift[k + n * XRES] = (1 / (double)XRES);
			dp_shift[k + n * XRES] = (1 / (double)XRES);
		}
	}


	_pdist_prior = nullptr;
	_pdist_lik = nullptr;
	_pdist_post = nullptr;

	_p_shift = nullptr;
	_dp_shift = nullptr;
	_dp = nullptr;
	_ddp = nullptr;

	if (!outputFile.is_open()) outputFile.open("output.csv", ios::out);
}

void TFcore::getSample() {
	samples++;

	// a vertical stack of row samples to compute MAV
	// the first row = the latest sample

	string data;
	string line;
	stringstream parse;

	//for (int n = 0; n < numChannels; n++) {
		//memcpy(&emgStack[1 + n * WIN_MAV], &emgStack[0 + n * WIN_MAV], (WIN_MAV - 1) * sizeof(double));
		//emgStack[n * WIN_MAV] = dataStack.front();
		//dataStack.erase(dataStack.begin());
	//}

	getline(inputFile, line, '\n');
	parse.str(line);
	for (int n = 0; n < numChannels; n++) {
		memcpy(&emgStack[1 + n * WIN_MAV], &emgStack[0 + n * WIN_MAV], (WIN_MAV - 1) * sizeof(double));
		getline(parse, data, ',');
		emgStack[n * WIN_MAV] = stod(data);
	}



	// compute MAV = a single row vector
	for (int n = 0; n < numChannels; n++) {
		emgMAV[n] = 0;
		for (int m = 0; m < WIN_MAV; m++) {
			emgMAV[n] += abs(emgStack[m + n * WIN_MAV]) / (double)(WIN_MAV);
		}
	}

	bCollect = true;
	bCompute = false;
}

void TFcore::proceedIteration() {
	if (isCollect() && (samples > 2 * WIN_MAV)) {
		bRegister = false;

		reglast += 1;

		p_max = -999;
		mpred = 0;

		for (int n = 0; n < numChannels; n++) {
			xnew[n] = (1 / xmax) * emgMAV[n];
			xnew[n] = min((double)1, xnew[n]);
			idnew[n] = floor(xnew[n] * (double)(XRES));
			if (idnew[n] > (XRES - 1)) idnew[n] = XRES - 1;
		}

		if (bFuncAdaptation) {
			if (bGPU_diffusion) {
				cudaMalloc(&_pdist_prior, M * XRES * numChannels * sizeof(double));
				cudaMalloc(&_p_shift, M * XRES * numChannels * sizeof(double));
				cudaMalloc(&_dp_shift, M * XRES * numChannels * sizeof(double));
				cudaMalloc(&_dp, M * XRES * numChannels * sizeof(double));
				cudaMalloc(&_ddp, M * XRES * numChannels * sizeof(double));

				// memcpy pdist_prior => pshift(dist_prior)
				for (int m = 0; m < M; m++) {
					for (int n = 0; n < numChannels; n++) {
						memcpy(&p_shift[1 + n * XRES + m * (XRES * numChannels)],
							&pdist_prior[0 + n * XRES + m * (XRES * numChannels)],
							(XRES - 1) * sizeof(double));
						p_shift[0 + n * XRES + m * (XRES * numChannels)] = p_shift[1 + n * XRES + m * (XRES * numChannels)];
					}
				}

				cudaMemcpy(_pdist_prior, pdist_prior, M * XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(_p_shift, p_shift, M * XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);

				_computeGradient << < 1, M* XRES* numChannels >> > (_pdist_prior, _p_shift, _dp, XRES);

				cudaMemcpy(dp, _dp, M * XRES * numChannels * sizeof(double), cudaMemcpyDeviceToHost);

				// memcpy dp(dist_prior) => dpshift(dist_prior)
				for (int m = 0; m < M; m++) {
					for (int n = 0; n < numChannels; n++) {
						memcpy(&dp_shift[1 + n * XRES + m * (XRES * numChannels)],
							&dp[0 + n * XRES + m * (XRES * numChannels)],
							(XRES - 1) * sizeof(double));
						dp_shift[0 + n * XRES + m * (XRES * numChannels)] = dp_shift[1 + n * XRES + m * (XRES * numChannels)];
					}
				}

				//cudaMemcpy(_dp, dp, M * XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(_dp_shift, dp_shift, M * XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);

				_computeGradient << < 1, M* XRES* numChannels >> > (_dp, _dp_shift, _ddp, XRES);

				cudaMemcpy(ddp, _ddp, M * XRES * numChannels * sizeof(double), cudaMemcpyDeviceToHost);

				_computeWeightedSum << < 1, M* XRES* numChannels >> > (_pdist_prior, _ddp, alpha);

				cudaMemcpy(pdist_prior, _pdist_prior, M * XRES * numChannels * sizeof(double), cudaMemcpyDeviceToHost);

				cudaFree(_pdist_prior);
				cudaFree(_p_shift);
				cudaFree(_dp_shift);
				cudaFree(_dp);
				cudaFree(_ddp);
			}
			else {
				computeDiffusion();
			}
			normalizePrior();

			computeLikelihood(xnew, sig2_update);
			//normalizeLikelihood();

			if (bGPU_product) {
				cudaMalloc(&_pdist_prior, M * XRES * numChannels * sizeof(double));
				cudaMalloc(&_pdist_lik, M * XRES * numChannels * sizeof(double));
				cudaMalloc(&_pdist_post, M * XRES * numChannels * sizeof(double));

				cudaMemcpy(_pdist_prior, pdist_prior, M * XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(_pdist_lik, pdist_lik, M * XRES * numChannels * sizeof(double), cudaMemcpyHostToDevice);

				_computeProduct << < M * XRES * numChannels, 1 >> > (_pdist_prior, _pdist_lik, _pdist_post);
				cudaMemcpy(pdist_post, _pdist_post, M * XRES * numChannels * sizeof(double), cudaMemcpyDeviceToHost);

				cudaFree(_pdist_prior);
				cudaFree(_pdist_lik);
				cudaFree(_pdist_post);
			}
			else {
				computeProduct();
			}
			normalizePost();

			// get maximum likelihood probability
			for (int m = 0; m < M; m++) {
				p_lik[m] = 0;
				for (int n = 0; n < numChannels; n++) {
					p_lik[m] += log10(pdist_post[idnew[n] + n * XRES + m * numChannels * XRES]);
				}
				//if (p_lik[m] < 1e-200) p_lik[m] = 1e-200;

				// MLE prediction
				if (p_lik[m] > p_max) {
					p_max = p_lik[m];
					mpred = m;
				}
			}
		}

		// registration
		if (bFuncRegistration) {
			if ((p_max < p_star) && (reglast > reghold) && (samples > 2 * WIN_MAV)) {
				registerPattern(xnew, sig2_reg);
				p_lik[M] = 0;
				for (int n = 0; n < numChannels; n++) {
					p_lik[M] += log10(pdist_post[idnew[n] + n * XRES + M * numChannels * XRES]);
				}
				p_max = p_lik[M];

				mpred = M;
				M += 1;
				reglast = 0;

				bRegister = true;
			}
			memcpy(&pdist_prior[mpred * XRES * numChannels], &pdist_post[mpred * XRES * numChannels], XRES * numChannels * sizeof(double));
		}

		bCollect = false;
		bCompute = true;
	}
	else {
		return;
	}
}

void TFcore::writeResult() {
	outputFile << p_max << ',' << mpred << ',';
	for (int n = 0; n < numChannels; n++) outputFile << xnew[n] << ',';

	outputFile << endl;
}

void TFcore::exportModel(fs::path filename) {
	fstream file;
	file.open(filename, ios::out | ios::binary);

	file.write(reinterpret_cast<char*>(&numChannels), sizeof(int));
	file.write(reinterpret_cast<char*>(&M_MAX), sizeof(int));
	file.write(reinterpret_cast<char*>(&xmax), sizeof(double));
	file.write(reinterpret_cast<char*>(&XRES), sizeof(int));
	file.write(reinterpret_cast<char*>(&WIN_MAV), sizeof(int));

	file.write(reinterpret_cast<char*>(&alpha), sizeof(double));
	file.write(reinterpret_cast<char*>(&beta), sizeof(double));
	file.write(reinterpret_cast<char*>(&p_star), sizeof(double));

	file.write(reinterpret_cast<char*>(&M), sizeof(int));
	file.write(reinterpret_cast<char*>(&reghold), sizeof(int));

	file.write(reinterpret_cast<char*>(&tableSize), sizeof(int));
	file.write(reinterpret_cast<char*>(&discretizeStep), sizeof(int));

	for (int k = 0; k < tableSize; k++) {
		file.write(reinterpret_cast<char*>(&normalTable[k]), sizeof(double));
	}

	for (int n = 0; n < numChannels; n++) {
		file.write(reinterpret_cast<char*>(&sig2_reg[n]), sizeof(double));
		file.write(reinterpret_cast<char*>(&sig2_update[n]), sizeof(double));
	}

	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			for (int k = 0; k < XRES; k++) {
				file.write(reinterpret_cast<char*>(&pdist_prior[k + n * XRES + m * XRES * numChannels]), sizeof(double));
				file.write(reinterpret_cast<char*>(&pdist_post[k + n * XRES + m * XRES * numChannels]), sizeof(double));
			}
		}
	}

	file.close();
}

void TFcore::importModel(fs::path filename) {
	fstream file;
	file.open(filename, ios::in | ios::binary);

	file.read(reinterpret_cast<char*>(&numChannels), sizeof(int));
	file.read(reinterpret_cast<char*>(&M_MAX), sizeof(int));
	file.read(reinterpret_cast<char*>(&xmax), sizeof(double));
	file.read(reinterpret_cast<char*>(&XRES), sizeof(int));
	file.read(reinterpret_cast<char*>(&WIN_MAV), sizeof(int));
	cout << "numChannels: " << numChannels << endl;
	cout << "M_MAX: " << M_MAX << endl;
	cout << "xmax: " << xmax << endl;
	cout << "XRES: " << XRES << endl;
	cout << "WIN_MAV: " << WIN_MAV << endl;

	file.read(reinterpret_cast<char*>(&alpha), sizeof(double));
	file.read(reinterpret_cast<char*>(&beta), sizeof(double));
	file.read(reinterpret_cast<char*>(&p_star), sizeof(double));
	cout << "alpha: " << alpha << endl;
	cout << "beta: " << beta << endl;
	cout << "p_star: " << p_star << endl;

	file.read(reinterpret_cast<char*>(&M), sizeof(int));
	file.read(reinterpret_cast<char*>(&reghold), sizeof(int));
	cout << "M: " << M << endl;
	cout << "reghold: " << reghold << endl;

	file.read(reinterpret_cast<char*>(&tableSize), sizeof(int));
	file.read(reinterpret_cast<char*>(&discretizeStep), sizeof(int));
	cout << "tableSize: " << tableSize << endl;
	cout << "discretizeStep: " << discretizeStep << endl;
	
	for (int k = 0; k < tableSize; k++) {
		file.read(reinterpret_cast<char*>(&normalTable[k]), sizeof(double));
	}
	cout << "normalTable" << endl;

	for (int n = 0; n < numChannels; n++) {
		file.read(reinterpret_cast<char*>(&sig2_reg[n]), sizeof(double));
		file.read(reinterpret_cast<char*>(&sig2_update[n]), sizeof(double));
	}
	cout << "sig2" << endl;

	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			for (int k = 0; k < XRES; k++) {
				file.read(reinterpret_cast<char*>(&pdist_prior[k + n * XRES + m * XRES * numChannels]), sizeof(double));
				file.read(reinterpret_cast<char*>(&pdist_post[k + n * XRES + m * XRES * numChannels]), sizeof(double));
			}
		}
	}
	cout << "pdist" << endl;

	file.close();
}

void TFcore::computeDiffusion() {
	//  pdist_prior => pshift(dist_prior)
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			memcpy(&p_shift[1 + n * XRES + m * (XRES * numChannels)],
				&pdist_prior[0 + n * XRES + m * (XRES * numChannels)],
				(XRES - 1) * sizeof(double));
			p_shift[0 + n * XRES + m * (XRES * numChannels)] = p_shift[1 + n * XRES + m * (XRES * numChannels)];
		}
	}

	// compute (p, pshift) => dp
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			for (int k = 0; k < XRES; k++) {
				dp[k + n * XRES + m * XRES * numChannels]
					= (double)XRES * (pdist_prior[k + n * XRES + m * XRES * numChannels] - p_shift[k + n * XRES + m * XRES * numChannels]);
			}
		}
	}

	// memcpy dp(dist_prior) => dpshift(dist_prior)
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			memcpy(&dp_shift[1 + n * XRES + m * (XRES * numChannels)],
				&dp[0 + n * XRES + m * (XRES * numChannels)],
				(XRES - 1) * sizeof(double));
			dp_shift[0 + n * XRES + m * (XRES * numChannels)] = dp_shift[1 + n * XRES + m * (XRES * numChannels)];
		}
	}

	// compute (dp, dpshift) => ddp
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < numChannels; n++) {
			for (int k = 0; k < XRES; k++) {
				ddp[k + n * XRES + m * XRES * numChannels]
					= (double)XRES * (dp[k + n * XRES + m * XRES * numChannels] - dp_shift[k + n * XRES + m * XRES * numChannels]);
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
		inputFile.open(filename.c_str(), ios::in);

		int m = 0;
		int n = 0;

		string data;
		string line;
		stringstream parse;

		getline(inputFile, line, '\n'); m++;
		parse.str(line);
		while (getline(parse, data, ',')) n++;

		while (getline(inputFile, line, '\n')) m++;
		inputFile.clear();
		inputFile.seekg(0, ios::beg);

		initModel(n);
		setNumSamples(m);
	}
}

void TFcore::loadData(fs::path filename) {
	if (filename.extension() == ".csv") {
		// comma separated file (samples x channels)
		inputFile.open(filename, ios::in);

		int m = 0;
		int n = 0;

		string data;
		string line;
		stringstream parse;

		getline(inputFile, line, '\n'); m++;
		parse.str(line);
		while (getline(parse, data, ',')) n++;

		while (getline(inputFile, line, '\n')) m++;
		inputFile.clear();
		inputFile.seekg(0, ios::beg);

		initModel(n);
		setNumSamples(m);
	}
}

void TFcore::loadData(vector<double> data, int ch) {
	initModel(ch);
	setNumSamples(data.size() / ch);
	for (int k = 0; k < data.size(); k++) dataStack.push_back(data[k]);
}

void TFcore::computeLikelihood(double* mu, double* sig2) {
	double* Z = (double*)calloc(XRES * numChannels, sizeof(double));
	int id = 0;
	for (int k = 0; k < XRES; k++) {
		for (int n = 0; n < numChannels; n++) {
			Z[k + n * XRES] = abs(((k + 1) / (double)XRES) - mu[n]);
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

	// memcpy 이전에 normalize 해버림
	double psum;
	for (int n = 0; n < numChannels; n++) {
		psum = 0;
		for (int k = 0; k < XRES; k++) psum += pdist_lik[k + n * XRES];
		for (int k = 0; k < XRES; k++) pdist_lik[k + n * XRES] /= psum;
	}

	if (M > 1) {
		for (int m = 1; m < M; m++) memcpy(&pdist_lik[m * XRES * numChannels], &pdist_lik[0], XRES * numChannels * sizeof(double));
	}
	free(Z);
}

void TFcore::registerPattern(double* mu, double* sig2) {
	double* Z = (double*)calloc(XRES * numChannels, sizeof(double));
	int id = 0;
	for (int n = 0; n < numChannels; n++) {
		for (int k = 0; k < XRES; k++) {
			Z[k + n * XRES] = abs(((k + 1) / (double)XRES) - mu[n]);
			Z[k + n * XRES] *= discretizeStep / sqrt(sig2[n]);

			id = floor(abs(Z[k + n * XRES]));
			if (id < tableSize) {
				pdist_post[k + n * XRES + M * XRES * numChannels] = normalTable[id];
			}
			else {
				pdist_post[k + n * XRES + M * XRES * numChannels] = EPSILON;
			}
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
