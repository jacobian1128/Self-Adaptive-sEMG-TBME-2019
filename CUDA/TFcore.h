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
	int M; // ���� �������� ��ü ������ ����
	int mpred; // classify�� ���� ������ �ε���

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
	double xmax; // ���͸��� sEMG ��ȣ�� �ִ��� ����

	Eigen::MatrixXd mu_init; // ù��° ������ ��� = �ַ� baseline noise amplitude
	Eigen::MatrixXd sig2_init; // ù��° ������ �л� = �ַ� baseline noise variance

	Eigen::MatrixXd sig2_reg; // ���ο� ������ ����� ���� �л�, �������� ���� ����
	Eigen::MatrixXd sig2_update; // ������ ���� �л�, �������� ���ϰ� �ݿ�

	double alpha; // Ȯ���� ����, Ŭ���� ���� Ȯ�� = ������ �а� ����
	double beta; // ������ �Ⱦ��� �Ķ����

	double pstar; // threshold probability, �� ���� ����� �ΰ���

	int regHold; // ��� ���Ŀ� ����� ���ø� ���� �� ����� ��������

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