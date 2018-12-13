close all
%% 1)Give sample time domain signal:
%Signal of ecg2
signal = load('ecg2.mat');
signal = signal.ecg2;

%Sample frequency and period
fs = 204.73;
Ts = 1/fs;

%Nyquist frequency
fn = fs / 2;

%Extrema and length of the raw signal
sig_max = max(signal);
sig_min = min(signal);
m = length(signal);
%Plot the signal in the time domain
totaltime = Ts*m;

time = linspace(0,totaltime,m);
plot(time,signal)
axis([0,totaltime,1.1*sig_min,1.1*sig_max]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2");

%% 2.a) Spectrum analysis:
figure

%Calculate the DFT of the signal using the FFT:
n=2^nextpow2(m); %Efficiently calculate FFT using N = power of 2 
X = fft(signal,n)/n; %MATLAB requires this get the correct amplitude
X_real = real(X); %Real component of the fft
X_imag = imag(X); %Imaginary component of the fft
X_abs = abs(X); %Absolute value of the fft

%Rescaling of the output:
%Rescaling based on reccomendations from https://nl.mathworks.com/help/matlab/ref/fft.html
f = fs*(0:(n/2))/n;
X_plot = X_abs(1:n/2+1);

%Plot the frequency spectrum
plot(f,X_plot)
title("Frequency spectrum of the ECG2");
legend("Absolute value");
xlabel("Frequentie (Hz)");

%% 2.b) Filter design
figure
%In the frequency spectrum, we powerline noise at 60Hz
%There appear to also be 2 harmonic compenents, at 180Hz(=3*60Hz) and at
%300Hz(=5*60Hz).
%We design 3 notch-filters to remove these 3 components.

%Define frequencies to be removed:
f0 = 50;
%only one notchfilter necessary because 3*50Hz=150Hz, this is more than the
%samplefrequency

%Angles on the unit circle for these frequencies:
theta0 = (pi/fn) * f0;

a = 0.9; %Value defined in lab assignment.

%Values for a-and b-component of the IIR. For the derivation of the
%formula, please refer to the lab report.
b0 = [1 -2*cos(theta0) 1];
a0 = [1 -2*a*cos(theta0) a*a];

%Placement of these poles and zeros on the unit circle:
zplane(b0,a0);

%% 3) Differential equations and direct form II schematic

%% 4) Impulse and frequency response
%These will only be calculated for the first filter (for 60Hz)
figure
freqz(b0,a0,1024);
fvtool(b0,a0)
%[H,W] = freqz(b0,a0,1024);
%plot(W/pi,20*log10(abs(H)));

figure
impz(b0,a0);
%% 5)Filtering:
signal0 = filter(b0,a0,signal);

%Spectrum-analysis after filtering:
figure
X0 = abs(fft(signal0,n) / n);
X0 = X0(1:n/2+1);

subplot(2,1,1)
plot(f,X_plot)
subplot(2,1,2)
plot(f,X0)

%Time-domain signals:
figure

subplot(2,1,1)
plot(time,signal)
time = linspace(0,totaltime,m);
axis([0,totaltime,1.1*sig_min,1.1*sig_max]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2");

subplot(2,1,2)
plot(time,signal0)
time = linspace(0,totaltime,m);
axis([0,totaltime,1.1*sig_min,1.1*sig_max]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2");

%% 5a) FIR Low-frequency drift removal
%Dit werkt al best goed(des te beter!): MA-filter met de IIR daarachter
figure
subplot(4,1,1)
plot(time,signal)
axis([0,totaltime,1.1*sig_min,1.1*sig_max]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2");

subplot(4,1,2)
b = [1/4, 1/4, 1/4, 1/4];
avgSignal = filter(b,1,signal);
plot(time,avgSignal)
axis([0,totaltime,1.1*min(avgSignal),1.1*max(avgSignal)]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2");

avgSignal = filter(b0,a0,avgSignal);

subplot(4,1,3)
plot(time,avgSignal)
axis([0,totaltime,1.1*sig_min,1.1*max(avgSignal)]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2");


%% BP Filter:
%We hebben nu 2 cutoff frequenties: wc1 en wc2: (als ge echt dringend moet)
f_pass1 = 40;
f_stop1 = 60;
f_pass2 = 0;
f_stop2 = 10;
w1s=f_stop1*pi;     %Stopband1: [w1s,pi]
w1p=f_pass1*pi;     %Passband1: [0,w1p]
w2p=f_pass2*pi;     %Passband2: [0,w2p]
w2s=f_stop2*pi;     %Stopband2: [w2s,pi]
% cut-off frequency halfway transition band
wc1=(w1p+w1s)/2;
wc2=(w2p+w2s)/2;

dw=min(w1p-w1s,w2s-w2p);


