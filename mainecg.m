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
subplot(1,3,1)
plot(absX)
title("FFT van signal1");
xlabel("Frequentie (Hz)");

%ingezoomde FFT
subplot(1,3,2);
axis([0,100,0,1.1*max(absX)]);
plot(absX(1:round(N/10)))
title("Ingezoomt tussen 0Hz en 100Hz");
xlabel("Frequentie (Hz)");
%==========================================================================
%Notch filter invoegen

%Notch filter zelf gemaakt
fs = 1000;
f0 = 60;
notchWidth = 0.9;
lnotch = notch(signal1, fs1, f0, notchWidth);
mnotch = notch(lnotch, fs1, f0, notchWidth); %geen idee waarom er twee keer wordt genotched hier
Xnotch = fft(mnotch);
absXnotch = abs(Xnotch);
hold on
%axis([0,100,0,1.1*max(abs(Xnotch))]);
plot(absXnotch(1:round(N/10)))
title("Ingezoomt tussen 0Hz en 100Hz, met notch filter");
xlabel("Frequentie (Hz)");
legend();

%FFT met notch maken
%subplot(1,3,3);
figure
axis([0,2500,1.1*ymin1,1.1*ymax1]); %waardes aanpassen
plot(absX)
title("FFT van signal1");
xlabel("Frequentie (Hz)");
hold on
plot(absXnotch)
legend();

%notch filter op signaal in de tijd
figure
plot(time1,signal1)
axis([0,totaltime1,1.1*ymin1,1.1*ymax1]);
title("Het originele signaal en het genotchte signaal");
xlabel("Tijd in s");
hold on
plot(time1, mnotch)
legend();

% %notch filter van matlab
% wo = 60/(fs1/2); %deze filter lijkt zelfs niet te werken 
% bw = wo/35;
% [c,d] = iirnotch(wo,bw);
% 
% nnotch = filter(c,d,signal1);
% figure 
% plot(time1,signal1)
% axis([0,totaltime1,1.1*ymin1,1.1*ymax1]);
% title("Het originele signaal en het 2de genotchte signaal");
% xlabel("Tijd in s");
% hold on
% plot(time1, nnotch)
% legend();
% 
% %FFT van matlabnotch volledig
% Ynotch = fft(nnotch);
% absYnotch = abs(Ynotch);
% figure
% plot(absX)
% title("FFT van signal1, met de notch van matlab");
% xlabel("Frequentie (Hz)");
% hold on
% plot(absYnotch)
% legend();
% 
% %FFT van matlabnotch ingezoomed
% figure
% axis([0,100,0,1.1*max(absX)]);
% plot(absX(1:round(N/10)))
% title("Ingezoomt tussen 0Hz en 100Hz");
% xlabel("Frequentie (Hz)");
% hold on
% %axis([0,100,0,1.1*max(abs(Xnotch))]);
% plot(absYnotch(1:round(N/10)))
% legend();

% %Notch filter zelf gedesigned
% %Ingezoomde FFT 
% onotch = designedNotch(signal1);
% Znotch = fft(onotch);
% absZnotch = abs(Znotch);
% figure
% axis([0,100,0,1.1*max(absX)]);
% plot(absX(1:round(N/10)))
% hold on
% plot(absZnotch(1:round(N/10)))
% title("Ingezoomt tussen 0Hz en 100Hz, met designed notch filter");
% xlabel("Frequentie (Hz)");
% legend();
% 
% %FFT met gedesignde notch maken
% figure
% plot(absX)
% title("FFT van signal1");
% xlabel("Frequentie (Hz)");
% hold on
% plot(absZnotch)
% legend();
% 
% %notch filter op signaal in de tijd
% figure
% plot(time1,signal1)
% axis([0,totaltime1,1.1*ymin1,1.1*ymax1]);
% title("Het originele signaal en het genotchte signaal");
% xlabel("Tijd in s");
% hold on
% plot(time1, onotch)
% legend();


