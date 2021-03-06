#include <cstdio>
#include <string>
#include <iostream>
#include <filesystem>

#include <thread>
#include <vector>

#include <Windows.h>

#include "TFcore.cuh"

using namespace std;

void worker(int iterSub, int iterCV) {
	char buff[100];
	sprintf(buff, "input/S%02d-train%d.csv", iterSub, iterCV);
	string fileInput(buff);
	fs::path pathInput(fileInput);

	sprintf(buff, "output0/S%02d-train%d.csv", iterSub, iterCV);
	string fileOutput(buff);
	fs::path pathOutput(fileOutput);

	//cout << fileInput << " >> " << fileOutput << endl;
	printf("S%02d CV%d :: %s >> %s\n", iterSub, iterCV, fileInput.c_str(), fileOutput.c_str());

	TFcore TF(pathOutput);
	TF.loadData(pathInput);
	//TF.resetModel();
	printf("S%02d CV%d :: data %d x %d\n", iterSub, iterCV, TF.getNumSamples(), TF.getNumChannels());

	//TF.enableGPU_Diffusion();
	//TF.enableGPU_Product();

	clock_t startTime = clock();
	while (TF.isRemaining()) {
		TF.getSample();
		TF.proceedIteration();
		TF.writeResult();

		if (TF.isRegister())
			printf("S%02d CV%d :: registration %d (%.2f %%)\n", iterSub, iterCV, TF.getNumPatterns(), 100 * (1 - TF.getNumRemaining() / (double)TF.getNumSamples()));
	}
	clock_t endTime = clock();
	sprintf(buff, "output0/S%02d-train%d-modelParam.bin", iterSub, iterCV);
	string modelOutput(buff);
	TF.exportModel(modelOutput);

	printf("S%02d CV%d :: %.2f [min]\n", iterSub, iterCV, (double)(endTime - startTime) / (60 * CLOCKS_PER_SEC));
}

wstring s2ws(const string& s)
{
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
	wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
	std::wstring r(buf);
	delete[] buf;
	return r;
}

void tgprintf(string hello) {
	string chat_id = "59233348";
	string token = "430183280:AAEoNtMwQQFxMhE7VP9hwQnHWI-RAtyD9No";
	string str = "https://api.telegram.org/bot" + token
		+ "/sendMessage?chat_id=" + chat_id + "&text=" + hello;
	
	wstring send = s2ws(str);
	printf("%s", str.c_str());
	ShellExecute(0, 0, send.c_str(), 0, SW_HIDE);
}


int main() {
	vector<thread> workers;
	tgprintf("Hello, World!");
	
	//int iterSub = 1;
	//int iterCV = 1;

	//for (iterSub = 4; iterSub < 5; iterSub++) {
	//	for (iterCV = 1; iterCV < 7; iterCV++) {
	//		workers.push_back(thread(worker, iterSub, iterCV));
	//	}
	//	for (int k = 0; k < workers.size(); k++) {
	//		workers[iterSub].join();
	//	}
	//	workers.clear();
	//}

	Sleep(10.00 * 1000);
	return 0;
}