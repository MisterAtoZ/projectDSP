close all
clc
% Eerste opdracht, de twee signalen plotten in de tijd
tabel = load('ecg.mat');
tabel2 = load('ecg.mat');

signal1 = tabel.ecg;
m1 = length(signal1);
ymax1 = max(signal1);
ymin1 = min(signal1);

signal2 = tabel2.ecg;
m2 = length(signal2);
ymax2 = max(signal2);
ymin2 = min(signal2);


fs1 = 1000;
Ts1 = 1/fs1;
fs2 = 204.73;
Ts2 = 1/fs2;

subplot(1,2,1)
totaltime1 = Ts1*m1;
totaltime2 = Ts2*m2;
time1 = linspace(0,totaltime1,m1);
plot(time1,signal1)
axis([0,totaltime1,1.1*ymin1,1.1*ymax1]);
xlabel("Time in s");

subplot(1,2,2)
time2 = linspace(0,totaltime2,m2);
plot(time2,signal2)
axis([0,totaltime2,1.1*ymin2,1.1*ymax2]);
xlabel("Time in s");

%FFT maken
X = fft(signal1);
Xreal = real(X);
Xcomp = imag(X);
absX = abs(X);
N = length(absX);
figure
subplot(1,2,1)
plot(absX)
title("FFT van signal1");
xlabel("Frequentie (Hz)");

%ingezoomde FFT
subplot(1,2,2);
axis([0,100,0,1.1*max(absX)]);
plot(absX(1:round(N/10)))
title("Ingezoomt tussen 0Hz en 100Hz");
xlabel("Frequentie (Hz)");


