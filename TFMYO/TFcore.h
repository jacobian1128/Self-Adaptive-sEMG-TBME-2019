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
	int M; // ���� �������� ��ü ������ ����
	int mpred; // classify�� ���� ������ �ε���

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
	double xmax; // ���͸��� sEMG ��ȣ�� �ִ��� ����

	MatrixXd mu_init; // ù��° ������ ��� = �ַ� baseline noise amplitude
	MatrixXd sig2_init; // ù��° ������ �л� = �ַ� baseline noise variance

	MatrixXd sig2_reg; // ���ο� ������ ����� ���� �л�, �������� ���� ����
	MatrixXd sig2_update; // ������ ���� �л�, �������� ���ϰ� �ݿ�

	double alpha; // Ȯ���� ����, Ŭ���� ���� Ȯ�� = ������ �а� ����
	double beta; // ������ �Ⱦ��� �Ķ����

	double pstar; // threshold probability, �� ���� ����� �ΰ���

	int regHold; // ��� ���Ŀ� ����� ���ø� ���� �� ����� ��������

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