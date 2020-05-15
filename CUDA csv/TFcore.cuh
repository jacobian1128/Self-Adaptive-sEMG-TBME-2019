#pragma once
#include <array>
#include <iostream>
#include <fstream>
#include <istream>

#define _USE_MATH_DEFINES
#include <math.h>
#include <iomanip>

#include <list>
#include <string>

#define WIN_MAV 256	// the number of MAV samples 
#define XRES 128	// resolution of discretizing probability distribution 
#define EPSILON 1.00e-10
#define M_MAX 2048

using namespace std;


class TFcore {
public:
	TFcore();
	~TFcore();

public:
	void initModel(int n);
	void loadData(string filename);

public:
	int getSamples() { return numsample; };
	int getPredict() { return mpred; };
	int getPatterns() { return M; };

	bool isCollect() { return bCollect; };
	bool isCompute() { return bCompute; };
	bool isRegister() { return bRegister; };
	bool isEmpty() { return (dataStack.size() == 0); };

private:
	int numch;
	int numsample;

	int reghold;
	double xmax;

	double alpha;
	double beta;

	double pstar;

	double* mu_init;
	double* sig2_init;

	double* sig2_reg;
	double* sig2_update;

private:
	bool bCollect;
	bool bCompute;
	bool bRegister;

private:
	double* pdist_prior;
	double* pdist_post;
	double* pdist_lik;

private:
	int mpred;
	int M;
	int reglast;
	double pmax;

	double* xbin;
	double* xnew;
	int* idnew;

	double* emgRaw;
	double* emgMAV;
	double* emgStack;

private:
	list<double> dataStack;

public:
	void getSample();
	void proceed();

private:
	void setSample(int k);

private:
	int tableSize;
	int discretizeStep;
	void computeNormal(double mu, double sig2, double* p);
	void computeNormal(double* mu, double* sig2, double* p);
	void normalizeProb(double* p);
	void computeDiffusion();

	void constructLookup();
	double* normalTable;

private:
	void setVal(double* vector, double val);
	void setZero(double* vector);
	void setOnes(double* vector);
};