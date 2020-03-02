#pragma once

#include <array>
#include <iostream>
#include <fstream>
#include <Windows.h>

#include <Eigen/Dense>

#define WIN_MAV 128	// the number of MAV samples 
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
	void getSample(std::array<int8_t, NUM_CH> emgSamples);
	void proceed();
	void initialize();

	void initDisplay();
	void printDisplay();

	bool isCollect();
	bool isCompute();
	bool isInitDisplay();
	
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

private:
	COORD posEMGText;
	COORD posParmText;
	COORD posProbText;

	void PrintText(const char * str, int x, int y, int tColor, int bColor);
	void PrintfXY(int x, int y, int tColor, int bColor, const char * str, ...);
	void PrintfXY(int x, int y, const char * str, ...);
	void printColorArea(int x, int y, int width, int height, int bColor);

};