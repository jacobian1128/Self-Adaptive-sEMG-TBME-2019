#ifndef TFcore_H_
#define TFcore_H_

#pragma once
#include <array>

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>

#define _USE_MATH_DEFINES
#include <math.h>
#include <iomanip>

#include <vector>
#include <string>

#define EPSILON 1.00e-10 // minimum value for probability distribution

using namespace std;
namespace fs = std::experimental::filesystem;

class CTFcoreMFCDlg;

class TFcore {
public:
	TFcore();
	TFcore(string filename);
	TFcore(fs::path filename);
	//TFcore(CTFcoreMFCDlg& dlg);
	~TFcore();

private:
	//CTFcoreMFCDlg& dlg_;

public:
	void loadData(string filename);
	void loadData(fs::path filename);
	void loadData(vector<double> data, int ch);

public:
	int getNumSamples() { return numSamples; }
	int getNumChannels() { return numChannels; }
	int getPredict() { return mpred; }
	int getNumPatterns() { return M; }
	//int getNumRemaining() { return ((int)dataStack.size() / numChannels); }
	int getNumRemaining() { return (numSamples - samples); }
	double getProbMax() { return p_max; }
	double getProbThresh() { return p_star; }

public:
	void resetModel() { initModel(numChannels); }
	void initModel(int n);
	void setNumSamples(int numsamples) { numSamples = numsamples; }
	void setNumChannels(int numchannels) { numChannels = numchannels; }

public:
	bool isCollect() { return bCollect; }
	bool isCompute() { return bCompute; }
	bool isRegister() { return bRegister; }

	bool isRemaining() { return !inputFile.eof(); }

	bool isGPU_Product() { return bGPU_product; }
	bool isGPU_Diffusion() { return bGPU_diffusion; }

	void enableGPU_Product() { bGPU_product = true; }
	void disableGPU_Product() { bGPU_product = false; }

	void enableGPU_Diffusion() { bGPU_diffusion = true; }
	void disableGPU_Diffusion() { bGPU_diffusion = false; }

	void enableRegistration() { bFuncRegistration = true; }
	void disableRegistration() { bFuncRegistration = false; }

	void enableAdaptation() { bFuncAdaptation = true; }
	void disableAdaptation() { bFuncAdaptation = false; }

public:
	void getSample();
	void proceedIteration();
	void writeResult();

public:
	void exportModel(fs::path filename);
	void importModel(fs::path filename);

private:
	int numChannels;
	int numSamples;

	int samples;
	int M_MAX;

	double xmax;

	int XRES;
	int WIN_MAV;

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
	bool bGPU_diffusion;

	bool bFuncRegistration;
	bool bFuncAdaptation;

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
	double* _pdist_prior;
	double* _pdist_lik;
	double* _pdist_post;

	double* _p_shift;
	double* _dp_shift;
	double* _dp;
	double* _ddp;

private:
	fstream inputFile;
	fstream outputFile;
};

#endif