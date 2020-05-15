#include "TFcore.h"

TFcore::TFcore() {
	// initialize model parameters = 손 댈일 있는 곳
	xmax = 2.00e-4;

	mu_init = 1.00e-2 * Eigen::MatrixXd::Ones(1, NUM_CH);
	sig2_init = 1.00e-3 * Eigen::MatrixXd::Ones(1, NUM_CH);

	sig2_reg = 1.00e-2 * Eigen::MatrixXd::Ones(1, NUM_CH);
	sig2_update = 1.00e+2 * Eigen::MatrixXd::Ones(1, NUM_CH);

	alpha = 1.00e-5;
	beta = 1.00e-10;

	pstar = 1e-17;


	regHold = 512;
	outputFile.open("data.txt");

	// initialize matrices = 손 댈일 없는 곳
	bCollect = false;
	bCompute = false;
	bInitDiaplay = false;
	bRegister = false;

	M = 1;
	regLast = 0;

	xbin = Eigen::MatrixXd::Zero(XRES, 1);
	for (int i = 0; i < XRES; i++) {
		xbin(i, 0) = (i / (double)(XRES));
	}
	emgRaw = Eigen::MatrixXd::Zero(1, NUM_CH);
	emgMAV = Eigen::MatrixXd::Zero(1, NUM_CH);
	emgStack = Eigen::MatrixXd::Zero(WIN_MAV, NUM_CH);

	xnew = Eigen::MatrixXd::Zero(1, NUM_CH);

	for (int m = 0; m < M_MAX; m++) {
		pdist_prior[m] = Eigen::MatrixXd::Zero(XRES, NUM_CH);
		pdist_post[m] = Eigen::MatrixXd::Zero(XRES, NUM_CH);
		mu[m] = Eigen::MatrixXd::Zero(1, NUM_CH);
	}
	pdist_lik = Eigen::MatrixXd::Zero(XRES, NUM_CH);
}

TFcore::~TFcore() {
	outputFile.close();
}

void TFcore::getSample(std::array<double, NUM_CH> emgSamples) {
	// a new sample = a row vector	
	for (int n = 0; n < NUM_CH; n++) {
		emgRaw(0, n) = abs((double)emgSamples[n]);
	}

	// a vertical stack of row samples to compute MAV
	emgStack.block(1, 0, WIN_MAV - 1, NUM_CH)
		= emgStack.block(0, 0, WIN_MAV - 1, NUM_CH);
	// the first row = the latest sample
	emgStack.block(0, 0, 1, NUM_CH) = emgRaw;

	// compute MAV = a single row vector
	emgMAV = emgStack.colwise().sum();
	emgMAV = (1 / (double)(WIN_MAV)) * emgMAV;

	bCollect = true;
	bCompute = false;
}

void TFcore::proceed() {
	bRegister = false;

	if (isCollect()) {
		regLast += 1;

		xnew = (1 / xmax) * emgMAV;
		for (int n = 0; n < NUM_CH; n++) {
			xnew(0, n) = min((double)1, xnew(0, n));
			idnew[n] = floor(xnew(0, n) * (double)(XRES));
			if (idnew[n] > (XRES - 1)) idnew[n] = XRES - 1;

			outputFile << xnew(0, n) << "\t";
		}
		outputFile << endl;

		computeDiffusion(); // pdist_prior = pdist_prior + alpha*dpdist_prior^2/dx^2
		pdist_lik = computeNormal(xnew, sig2_update);

		pmax = -1;
		mpred = 0;

		// 모든 패턴들에 대해 평가함
		for (int m = 0; m < M; m++) {
			//pdist_post[m] = pdist_prior[m];
			pdist_post[m] = pdist_prior[m].array() * pdist_lik.array();
			pdist_post[m] = normalizeProb(pdist_post[m]);

			// 각 패턴별 likelihood probability를 계산함
			plik[m] = 1;
			for (int n = 0; n < NUM_CH; n++) {
				plik[m] *= pdist_post[m](idnew[n], n);
			}
			if (plik[m] < 1e-99) plik[m] = 1e-99;

			if (plik[m] > pmax) {
				pmax = plik[m];
				mpred = m;
			}
		}

		// 새 패턴 등록하기
		if ((pmax < pstar) && (regLast > regHold)) {
			pdist_prior[M] = computeNormal(xnew, sig2_reg);
			pdist_post[M] = pdist_prior[M];
			plik[M] = 1;
			for (int n = 0; n < NUM_CH; n++) {
				plik[M] *= pdist_post[M](idnew[n], n);
			}

			mpred = M;
			M += 1;
			regLast = 0;

			bRegister = true;
		}

		for (int m = 0; m < M; m++) {
			//pdist_prior[m] = pdist_post[m];
			mu[m] = computeMean(pdist_prior[m]);
		}

		pdist_prior[mpred] = pdist_post[mpred];
		//pdist_prior[mpred] = pdist_post[mpred].array() * pdist_lik.array();

		bCollect = false;
		bCompute = true;
	}
}

void TFcore::initialize() {
	int m = 0;
	pdist_prior[m] = computeNormal(mu_init, sig2_init);
	mu[m] = computeMean(pdist_prior[m]);

	constructLookup();
}








Eigen::MatrixXd TFcore::computeNormal(double mu, double sig2) {
	Eigen::MatrixXd prob(XRES, 1);
	// 	prob = xbin - mu * MatrixXd::Ones(XRES, 1);
	// 	prob = (-1 / (2 * sig2)) * prob.array().square();
	// 	prob = (1 / sqrt(2 * M_PI * sig2)) * prob.array().exp();

	Eigen::MatrixXd Z(XRES, 1);
	Z = xbin - mu * Eigen::MatrixXd::Ones(XRES, 1);
	Z = 100 * Z * (1 / sqrt(sig2));
	for (int i = 0; i < XRES; i++) {
		if (floor(abs(Z(i, 0))) < normalTable.cols())
			prob(i, 0) = normalTable(1, floor(abs(Z(i, 0))));
		else
			prob(i, 0) = 0;
	}

	prob = normalizeProb(prob);
	return prob;
}

Eigen::MatrixXd TFcore::computeNormal(Eigen::MatrixXd mu, Eigen::MatrixXd sig2) {
	Eigen::MatrixXd prob(XRES, NUM_CH);
	for (int n = 0; n < NUM_CH; n++) {
		prob.block(0, n, XRES, 1) = computeNormal(mu(0, n), sig2(0, n));
	}
	return prob;
}

Eigen::MatrixXd TFcore::computeMean(Eigen::MatrixXd prob) {
	Eigen::MatrixXd mu(1, NUM_CH);
	mu = xbin.transpose() * prob;
	return mu;
}

Eigen::MatrixXd TFcore::normalizeProb(Eigen::MatrixXd prob) {
	for (int n = 0; n < prob.cols(); n++) {
		prob.block(0, n, XRES, 1) = (1 / prob.block(0, n, XRES, 1).array().sum()) * prob.block(0, n, XRES, 1);
	}
	return prob;
}

void TFcore::computeDiffusion() {
	Eigen::MatrixXd df = Eigen::MatrixXd::Zero(XRES, NUM_CH);
	Eigen::MatrixXd ddf = Eigen::MatrixXd::Zero(XRES, NUM_CH);
	for (int m = 0; m < M; m++) {
		df.block(0, 0, XRES - 1, NUM_CH) = pdist_prior[m].block(1, 0, XRES - 1, NUM_CH) - pdist_prior[m].block(0, 0, XRES - 1, NUM_CH);
		df.block(XRES - 1, 0, 1, NUM_CH) = df.block(XRES - 2, 0, 1, NUM_CH);
		ddf.block(0, 0, XRES - 1, NUM_CH) = df.block(1, 0, XRES - 1, NUM_CH) - df.block(0, 0, XRES - 1, NUM_CH);
		ddf.block(XRES - 1, 0, 1, NUM_CH) = ddf.block(XRES - 2, 0, 1, NUM_CH);
		pdist_prior[m] = pdist_prior[m] + (alpha / (XRES * XRES)) * ddf;
	}
}

void TFcore::constructLookup() {
	normalTable = Eigen::MatrixXd::Zero(1, 500);
	for (int i = 0; i < normalTable.cols(); i++) {
		normalTable(0, i) = (1 / sqrt(2 * M_PI)) * exp(-0.5 * ((double)i / 100) * ((double)i / 100));
	}
}