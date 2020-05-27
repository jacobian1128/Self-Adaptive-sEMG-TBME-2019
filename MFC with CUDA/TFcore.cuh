#ifndef TFcore_H_
#define TFcore_H_

#pragma once
#include <array>
#include <iostream>
#include <fstream>
#include <istream>

#define _USE_MATH_DEFINES
#include <math.h>
#include <iomanip>

#include <vector>
#include <string>

#define WIN_MAV 2	// the number of MAV samples 
#define XRES 128	// resolution of discretizing probability distribution 
#define EPSILON 1.00e-10 // minimum value for probability distribution
#define M_MAX 2048

using namespace std;

class CTFcoreMFCDlg;

class TFcore {
public:
	TFcore(CTFcoreMFCDlg& dlg);
	~TFcore();

private:
	CTFcoreMFCDlg& dlg_;

public:
	void loadData(string filename);
	void loadData(vector<double> data, int ch);

public:
	int getNumSamples() { return numSamples; }
	int getNumChannels() { return numChannels; }
	int getPredict() { return mpred; }
	int getNumPatterns() { return M; }
	int getRemaining() { return ((int)dataStack.size() / numChannels); }
	double getProbMax() { return p_max; }
	double getProbThresh() { return p_star; }

public:
	void initModel(int n);
	void setNumSamples(int numsamples) { numSamples = numsamples; }
	void setNumChannels(int numchannels) { numChannels = numchannels; }

public:
	bool isCollect() { return bCollect; }
	bool isCompute() { return bCompute; }
	bool isRegister() { return bRegister; }

	bool isGPU_Product() { return bGPU_product; }

	void enableGPU_Product() { bGPU_product = true; }
	void disableGPU_Product() { bGPU_product = false; }

public:
	void getSample();
	void proceedIteration();
	void writeResult();

private:
	int numChannels;
	int numSamples;

	double xmax;

	double alpha;
	double beta;

	double p_star;
	double* p_lik;

	double* mu_init;
	double* sig2_init;

	double* sig2_reg;
	double* sig2_update;

private:
	bool bCollect;
	bool bCompute;
	bool bRegister;

private:
	bool bGPU_product;

public:
	double* pdist_prior;
	double* pdist_post;
	double* pdist_lik;

	double* dp;
	double* ddp;
	double* p_shift;
	double* dp_shift;

private:
	int mpred;
	int M;
	int reglast;
	int reghold;
	double p_max;

	double* xbin;
	double* xnew;
	int* idnew;

	double* emgMAV;
	double* emgStack;

private:
	vector<double> dataStack;

private:
	int tableSize;
	int discretizeStep;

	void registerPattern(double* mu, double* sig2);
	void computeLikelihood(double* mu, double* sig2);
	void computeProduct();
	void computeDiffusion();

	void normalizeLikelihood();
	void normalizePrior();
	void normalizePost();

	double* normalTable;
	void constructLookup();

private:
	void setVal(double* vector, double val);
	void setZero(double* vector);
	void setOnes(double* vector);

private:
	fstream outFile;
};

#endif