#include <array>
#include <iostream>
#include <fstream>
#include <Windows.h>

#include <math.h>

#include "Eigen/Dense"

#include "mex.hpp"
#include "mexAdapter.hpp"

#define _USE_MATH_DEFINES

#define WIN_MAV 32	// the number of MAV samples 
#define NUM_CH 8	// the number of sEMG channels
#define XRES 128	// resolution of discretizing probability distribution 
#define M_MAX 128	// maximum number of patterns

using namespace Eigen;
using namespace std;

class TFcore {
public:
	TFcore();
	~TFcore();

public:
	void getSample(std::array<double, NUM_CH> emgSamples);
	void proceed();
	void initialize();

	bool isCollect();
	bool isCompute();

	int getPredict() { return mpred; };

private:
	bool bCollect;
	bool bCompute;
	bool bInitDiaplay;

private:
	MatrixXd emgRaw;
	MatrixXd emgMAV;
	MatrixXd emgStack;

	MatrixXd xnew;

private:
	ofstream outputFile;

private:
	int M; // 현재 시점에서 전체 패턴의 갯수
	int mpred; // classify된 현재 패턴의 인덱스

	int idnew[NUM_CH];
	double pmax;
	int regLast;

private:
	MatrixXd xbin;
	MatrixXd pdist_prior[M_MAX];
	MatrixXd pdist_post[M_MAX];
	MatrixXd pdist_lik;

	MatrixXd mu[M_MAX];
	double plik[M_MAX];

private:
	double xmax; // 필터링된 sEMG 신호의 최댓값을 설정

	MatrixXd mu_init; // 첫번째 패턴의 평균 = 주로 baseline noise amplitude
	MatrixXd sig2_init; // 첫번째 패턴의 분산 = 주로 baseline noise variance

	MatrixXd sig2_reg; // 새로운 패턴을 등록할 때의 분산, 낮을수록 좁은 분포
	MatrixXd sig2_update; // 갱신할 때의 분산, 낮을수록 강하게 반영

	double alpha; // 확산의 정도, 클수록 강한 확산 = 분포가 넓게 퍼짐
	double beta; // 지금은 안쓰는 파라미터

	double pstar; // threshold probability, 새 패턴 등록의 민감도

	int regHold; // 등록 이후에 몇번의 샘플링 동안 새 등록을 금지할지

private:
	MatrixXd computeNormal(double mu, double sig2);
	MatrixXd computeNormal(MatrixXd mu, MatrixXd sig2);
	MatrixXd computeMean(MatrixXd prob);
	MatrixXd normalizeProb(MatrixXd prob);
	void computeDiffusion();
};


class MexFunction : public matlab::mex::Function {
    
public:
	TFcore TF;
	std::array<double, NUM_CH> emgSamples;

	std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    std::ostringstream stream;
    matlab::data::ArrayFactory factory;
    
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        matlab::data::TypedArray<double> inArray = inputs[0];
        
		matlab::data::ArrayDimensions dim = inArray.getDimensions();
        stream << "ArrayDimensions: " << dim[0] << " " << dim[1] << std::endl;
		flushStream();

        dim[1] = 1;
        matlab::data::TypedArray<int32_t> outArray = factory.createArray<int32_t>(dim);
        
		for (int i = 0; i < dim[0]; i++) {
            for (int ch = 0; ch < NUM_CH; ch++) {
                emgSamples[ch] = inArray[i][ch];
            }
			TF.getSample(emgSamples);
			TF.proceed();

			outArray[i] = (int32_t)TF.getPredict();
        }
        outputs[0] = outArray;
        
    }
    
public:
    void flushStream() {
        matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({ factory.createScalar(stream.str()) }));
        stream.str("");
    }
};


TFcore::TFcore() {
	// initialize model parameters = 손 댈일 있는 곳
	xmax = 1.00;

	mu_init = 1.00e-2 * MatrixXd::Ones(1, NUM_CH);
	sig2_init = 1.00e-3 * MatrixXd::Ones(1, NUM_CH);

	sig2_reg = 1.00e-2 * MatrixXd::Ones(1, NUM_CH);
	sig2_update = 1.00e+2 * MatrixXd::Ones(1, NUM_CH);

	alpha = 1.00e-5;
	beta = 1.00e-10;

	pstar = 1e-24;


	regHold = 128;
	outputFile.open("data.txt");

	// initialize matrices = 손 댈일 없는 곳
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
		// 모든 패턴들에 대해 평가함
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
	int m = 0;
	pdist_prior[m] = computeNormal(mu_init, sig2_init);
	mu[m] = computeMean(pdist_prior[m]);
}

bool TFcore::isCollect() { return bCollect; }
bool TFcore::isCompute() { return bCompute; }








MatrixXd TFcore::computeNormal(double mu, double sig2) {
	MatrixXd prob(XRES, 1);
	prob = xbin - mu * MatrixXd::Ones(XRES, 1);
	prob = (-1 / (2 * sig2)) * prob.array().square();
	prob = (1 / sqrt(2 * M_PI * sig2)) * prob.array().exp();

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
		pdist_prior[m] = pdist_prior[m] + (alpha / (XRES * XRES)) * ddf;
	}
}