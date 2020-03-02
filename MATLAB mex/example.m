clear;

load('data.mat');
mex computeTFMEX.cpp;

predict = computeTFMEX(emg);

figure(1);
subplot(2,1,1); plot(time,emg);
subplot(2,1,2); plot(time,predict);
