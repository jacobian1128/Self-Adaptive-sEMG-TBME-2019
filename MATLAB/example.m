clear;

load('data.mat');
mex computeTFMEX.cpp;

predict = computeTFMEX(emg);

figure(1);
subplot(3,1,1); plot(time,gyro);
xlabel('time [s]'); ylabel('angular velocity [deg/s]');
subplot(3,1,2); plot(time,emg);
xlabel('time [s]'); ylabel('amplified sEMG [V]');
subplot(3,1,3); plot(time,predict);
xlabel('time [s]'); ylabel('classified pattern');
print('figure.png','-dpng');