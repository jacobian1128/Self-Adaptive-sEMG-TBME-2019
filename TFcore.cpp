#define _USE_MATH_DEFINES

#include <math.h>

#include "WinConsoleCtrl.h"
#include "TFcore.h"

TFcore::TFcore() {
	// initialize model parameters
	xmax = 80;

	mu_init = 1.00e-2*MatrixXd::Ones(1, NUM_CH);
	sig2_init = 1.00e-3*MatrixXd::Ones(1, NUM_CH);

	sig2_reg = 1.00e-2*MatrixXd::Ones(1, NUM_CH);
	sig2_update = 1.00e+2*MatrixXd::Ones(1, NUM_CH);

	alpha = 1.00e-5;
	beta = 1.00e-10;

	pstar = 1e-24;


	regHold = 128;
	outputFile.open("data.txt");
	
	// initialize matrices (do not modify)
	bCollect = false;
	bCompute = false;
	bInitDiaplay = false;

	M = 1;
	regLast = 0;

	xbin = MatrixXd::Zero(XRES, 1);
	for (int i = 0; i < XRES; i++) {
		xbin(i, 0) = (i / (double)(XRES));
	}
	emgRaw = MatrixXd::Zero(1, NUM_CH);
	emgMAV = MatrixXd::Zero(1, NUM_CH);
	emgStack = MatrixXd::Zero(WIN_MAV, NUM_CH);

	xnew = MatrixXd::Zero(1, NUM_CH);

	for (int m = 0; m < M_MAX; m++) {
		pdist_prior[m] = MatrixXd::Zero(XRES, NUM_CH);
		pdist_post[m] = MatrixXd::Zero(XRES, NUM_CH);
		mu[m] = MatrixXd::Zero(1, NUM_CH);
	}
	pdist_lik = MatrixXd::Zero(XRES, NUM_CH);
}

TFcore::~TFcore() {
	outputFile.close();
}

void TFcore::getSample(std::array<int8_t, NUM_CH> emgSamples) {
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
	if (isCollect()) {
		regLast += 1;

		xnew = (1 / xmax) * emgMAV;
		for (int n = 0; n < NUM_CH; n++) {
			xnew(0, n) = min(1, xnew(0, n));
			idnew[n] = floor(xnew(0, n) * (double)(XRES));
			if (idnew[n] > (XRES - 1)) idnew[n] = XRES - 1;

			outputFile << xnew(0, n) << "\t";
		}
		outputFile << endl;

		computeDiffusion(); // pdist_prior = pdist_prior + alpha*dpdist_prior^2/dx^2
		pdist_lik = computeNormal(xnew, sig2_update);

		pmax = -1;
		mpred = 0;
		// evaluate for all distributions
		for (int m = 0; m < M; m++) {
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

		// register a new pattern
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

			initDisplay();
		}

		for (int m = 0; m < M; m++) {
			//pdist_prior[m] = pdist_post[m];
			mu[m] = computeMean(pdist_prior[m]);
		}

		pdist_prior[mpred] = pdist_post[mpred];

		bCollect = false;
		bCompute = true;
	}
}

void TFcore::initialize() {
	system("cls");

	int m = 0;
	pdist_prior[m] = computeNormal(mu_init, sig2_init);
	mu[m] = computeMean(pdist_prior[m]);
}

bool TFcore::isCollect() { return bCollect; }
bool TFcore::isCompute() { return bCompute; }
bool TFcore::isInitDisplay() { return bInitDiaplay; }








MatrixXd TFcore::computeNormal(double mu, double sig2) {
	MatrixXd prob(XRES, 1);
	prob = xbin - mu*MatrixXd::Ones(XRES, 1);
	prob = (-1 / (2*sig2)) * prob.array().square();
	prob = (1 / sqrt(2*M_PI*sig2)) * prob.array().exp();

	prob = normalizeProb(prob);
	return prob;
}

MatrixXd TFcore::computeNormal(MatrixXd mu, MatrixXd sig2) {
	MatrixXd prob(XRES, NUM_CH);
	for (int n = 0; n < NUM_CH; n++) {
		prob.block(0, n, XRES, 1) = computeNormal(mu(0, n), sig2(0, n));
	}
	return prob;
}

MatrixXd TFcore::computeMean(MatrixXd prob) {
	MatrixXd mu(1, NUM_CH);
	mu = xbin.transpose() * prob;
	return mu;
}

MatrixXd TFcore::normalizeProb(MatrixXd prob) {
	for (int n = 0; n < prob.cols(); n++) {
		prob.block(0, n, XRES, 1) = (1 / prob.block(0, n, XRES, 1).array().sum()) * prob.block(0, n, XRES, 1);
	}
	return prob;
}

void TFcore::computeDiffusion() {
	MatrixXd df = MatrixXd::Zero(XRES, NUM_CH);
	MatrixXd ddf = MatrixXd::Zero(XRES, NUM_CH);
	for (int m = 0; m < M; m++) {
		df.block(0, 0, XRES - 1, NUM_CH) = pdist_prior[m].block(1, 0, XRES - 1, NUM_CH) - pdist_prior[m].block(0, 0, XRES - 1, NUM_CH);
		df.block(XRES - 1, 0, 1, NUM_CH) = df.block(XRES - 2, 0, 1, NUM_CH);
		ddf.block(0, 0, XRES - 1, NUM_CH) = df.block(1, 0, XRES - 1, NUM_CH) - df.block(0, 0, XRES - 1, NUM_CH);
		ddf.block(XRES - 1, 0, 1, NUM_CH) = ddf.block(XRES - 2, 0, 1, NUM_CH);
		pdist_prior[m] = pdist_prior[m] + (alpha / (XRES * XRES))*ddf;
	}
}







void TFcore::initDisplay() {
	bInitDiaplay = false;
	ShowCCursor(false);
	SetConsoleSize(120, 50, 120, 30);

	system("cls");
	PrintfXY(0, 0, "");

	//	   "01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789"
	//	   "0         1         2         3         4         5         6         7         8         9         A"
	posParmText = GetCCursorPos(); posParmText.Y += 2;
	printf("##################################################################################################\n");
	printf("                            model parameter                                                       \n");
	printf("    xmax 00   sig2 init 0.00E+00  reg 0.00E+00  update 0.00E+00   alpha 0.00E+00 beta 0.00E+00    \n");

	posEMGText = GetCCursorPos(); posEMGText.Y += 2;
	printf("##################################################################################################\n");
	printf(" index   | Sen1       Sen2       Sen3       Sen4       Sen5       Sen6       Sen7       Sen8      \n");
	printf(" raw     |   000.00     000.00     000.00     000.00     000.00     000.00     000.00     000.00  \n");
	printf(" norm    |   0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000  \n");
	//printf(" idnew   |      000        000        000        000        000        000        000        000  \n");

	posProbText = GetCCursorPos(); posProbText.Y += 1;
	printf("##################################################################################################\n");
	for (int m = 0; m < M; m++) {
	printf(" ptrn %2d |   0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000   0.00E+00\n", m);
	}
	bInitDiaplay = true;
}

void TFcore::printDisplay() {
	if (isInitDisplay()) {
		int x, y;

		// display model parameters
		x = 0;
		y = posParmText.Y;
		PrintfXY(x + 9, y, "%2.0f", xmax);
		PrintfXY(x + 24, y, "%8.2E", sig2_init(0, 0));
		PrintfXY(x + 38, y, "%8.2E", sig2_reg(0, 0));
		PrintfXY(x + 55, y, "%8.2E", sig2_update(0, 0));
		PrintfXY(x + 72, y, "%8.2E", alpha);
		PrintfXY(x + 86, y, "%8.2E", beta);

		y = posProbText.Y - 1;
		PrintfXY(x + 99, y, "%8.2E", pstar);

		// display sEMG values
		for (int n = 0; n < NUM_CH; n++) {
			x = 11 * (n + 1);
			y = posEMGText.Y;
			PrintfXY(x, y + 0, "  %6.2f", emgMAV(0, n));

			if (xnew(0, n) < 1) PrintfXY(x, y + 1, "  %6.4f", xnew(0, n));
			else				PrintfXY(x, y + 1, CC_RED, CC_DEFAULT, "  %6.4f", xnew(0, n));

			//PrintfXY(x + 3, y + 2, "  %3d", idnew[n]);
		}

		// display patterns
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < NUM_CH; n++) {
				x = 11 * (n + 1);
				y = posProbText.Y + m;
				PrintfXY(x, y, "  %6.4f", mu[m](0, n));
			}
			PrintfXY(0, y, " ptrn %02d ", m);
			PrintfXY(99, y, "%8.2E", plik[m]);
		}

		PrintfXY(0, posProbText.Y + mpred, CC_BLACK, CC_RED, " ptrn %02d ", mpred);
		PrintfXY(99, posProbText.Y + mpred, CC_BLACK, CC_RED, "%8.2E", plik[mpred]);
	}
}

void TFcore::PrintText(const char * str, int x, int y, int tColor, int bColor)
{
	int prevTextColor = 0;
	int prevBkgndColor = 0;

	GetCColor(&prevTextColor, &prevBkgndColor);

	SetCCursorPos(x, y);
	SetCColor(tColor, bColor);

	printf(str);

	SetCColor(prevTextColor, prevBkgndColor);
}

void TFcore::PrintfXY(int x, int y, int tColor, int bColor, const char * str, ...)
{
	va_list args;

	int prevTextColor = 0;
	int prevBkgndColor = 0;
	GetCColor(&prevTextColor, &prevBkgndColor);

	SetCCursorPos(x, y);
	SetCColor(tColor, bColor);

	va_start(args, str);
	vprintf(str, args);
	va_end(args);

	SetCColor(prevTextColor, prevBkgndColor);
}

void TFcore::PrintfXY(int x, int y, const char * str, ...)
{
	va_list args;

	SetCCursorPos(x, y);

	va_start(args, str);
	vprintf(str, args);
	va_end(args);
}

void TFcore::printColorArea(int x, int y, int width, int height, int bColor)
{
	int prevBkgndColor = GetCBkgndColor();
	SetCCursorPos(x, y);
	SetCBkgndColor(bColor);
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++) printf(" ");
		printf("\n");
	}
	SetCBkgndColor(prevBkgndColor);
}
