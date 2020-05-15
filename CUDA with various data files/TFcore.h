#pragma once

#include <array>
#include <iostream>
#include <fstream>

#include <math.h>
#include <iomanip>

#include "Eigen/Dense"

#define _USE_MATH_DEFINES

#define WIN_MAV 256	// the number of MAV samples 
#define NUM_CH 12	// the number of sEMG channels
#define XRES 128	// resolution of discretizing probability distribution 
#define M_MAX 128	// maximum number of patterns

// using namespace Eigen;
// using namespace std;

class TFcore {
public:
	TFcore();
	~TFcore();

public:
	void getSample(std::array<double, NUM_CH> emgSamples);
	void proceed();
	void initialize();

	bool isCollect() { return bCollect; };
	bool isCompute() { return bCompute; };
	bool isRegister() { return bRegister; };

	int getPredict() { return mpred; };
	int getPatterns() { return M; };

private:
	bool bCollect;
	bool bCompute;
	bool bInitDiaplay;
	bool bRegister;

private:
	Eigen::MatrixXd emgRaw;
	Eigen::MatrixXd emgMAV;
	Eigen::MatrixXd emgStack;

	Eigen::MatrixXd xnew;

private:
	std::ofstream outputFile;

private:
	int M; // 현재 시점에서 전체 패턴의 갯수
	int mpred; // classify된 현재 패턴의 인덱스

	int idnew[NUM_CH];
	double pmax;
	int regLast;

private:
	Eigen::MatrixXd xbin;
	Eigen::MatrixXd pdist_prior[M_MAX];
	Eigen::MatrixXd pdist_post[M_MAX];
	Eigen::MatrixXd pdist_lik;

	Eigen::MatrixXd mu[M_MAX];
	double plik[M_MAX];

private:
	double xmax; // 필터링된 sEMG 신호의 최댓값을 설정

	Eigen::MatrixXd mu_init; // 첫번째 패턴의 평균 = 주로 baseline noise amplitude
	Eigen::MatrixXd sig2_init; // 첫번째 패턴의 분산 = 주로 baseline noise variance

	Eigen::MatrixXd sig2_reg; // 새로운 패턴을 등록할 때의 분산, 낮을수록 좁은 분포
	Eigen::MatrixXd sig2_update; // 갱신할 때의 분산, 낮을수록 강하게 반영

	double alpha; // 확산의 정도, 클수록 강한 확산 = 분포가 넓게 퍼짐
	double beta; // 지금은 안쓰는 파라미터

	double pstar; // threshold probability, 새 패턴 등록의 민감도

	int regHold; // 등록 이후에 몇번의 샘플링 동안 새 등록을 금지할지

private:
	Eigen::MatrixXd computeNormal(double mu, double sig2);
	Eigen::MatrixXd computeNormal(Eigen::MatrixXd mu, Eigen::MatrixXd sig2);
	Eigen::MatrixXd computeMean(Eigen::MatrixXd prob);
	Eigen::MatrixXd normalizeProb(Eigen::MatrixXd prob);
	void computeDiffusion();

private:
	void constructLookup();
	Eigen::MatrixXd normalTable;
};