close all
%% 1)Give sample time domain signal:
%Signal of ecg1
signal = load('ecg.mat');
signal = signal.ecg;

%Sample frequency and period
fs = 1000;
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
title("ECG");

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
title("Frequency spectrum of the ECG");
legend("Absolute value");
xlabel("Frequentie (Hz)");

%% 2.b) Filter design
figure
%In the frequency spectrum, we powerline noise at 60Hz
%There appear to also be 2 harmonic compenents, at 180Hz(=3*60Hz) and at
%300Hz(=5*60Hz).
%We design 3 notch-filters to remove these 3 components.

%Define frequencies to be removed:
f0 = 60;
f1 = 3*60;
f2 = 5*60;

%Angles on the unit circle for these frequencies:
theta0 = (pi/fn) * f0;
theta1 = (pi/fn) * f1;
theta2 = (pi/fn) * f2;

a = 0.9; %Value defined in lab assignment.

%Values for a-and b-component of the IIR. For the derivation of the
%formula, please refer to the lab report.
b0 = [1 -2*cos(theta0) 1];
b1 = [1 -2*cos(theta1) 1];
b2 = [1 -2*cos(theta2) 1];

a0 = [1 -2*a*cos(theta0) a*a];
a1 = [1 -2*a*cos(theta1) a*a];
a2 = [1 -2*a*cos(theta2) a*a];

%Placement of these poles and zeros on the unit circle:
subplot(1,3,1);
zplane(b0,a0);
subplot(1,3,2);
zplane(b1,a1);
subplot(1,3,3);
zplane(b2,a2);

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
signal1 = filter(b1,a1,signal0);
signal2 = filter(b2,a2,signal1);

%Spectrum-analysis after filtering:
figure
X0 = abs(fft(signal0,n) / n);
X1 = abs(fft(signal1,n) / n);
X2 = abs(fft(signal2,n) / n);

X0 = X0(1:n/2+1);
X1 = X1(1:n/2+1);
X2 = X2(1:n/2+1);
subplot(4,1,1)
plot(f,X_plot)
xlabel("Frequency in Hertz");
title("original ECG");
subplot(4,1,2)
plot(f,X0)
xlabel("Frequency in Hertz");
title("60Hz notch filter ECG");
subplot(4,1,3)
plot(f,X1)
xlabel("Frequency in Hertz");
title("60Hz and 180Hz notch filter ECG");
subplot(4,1,4)
plot(f,X2)
xlabel("Frequency in Hertz");
title("60Hz, 180Hz and 300Hz notch filter ECG");

%Time-domain signals:
figure

subplot(4,1,1)
plot(time,signal)
time = linspace(0,totaltime,m);
axis([0,totaltime,1.1*sig_min,1.1*sig_max]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("Originele ECG");

subplot(4,1,2)
plot(time,signal0)
time = linspace(0,totaltime,m);
axis([0,totaltime,1.1*sig_min,1.1*sig_max]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("60Hz notch filter ECG");

subplot(4,1,3)
plot(time,signal1)
time = linspace(0,totaltime,m);
axis([0,totaltime,1.1*sig_min,1.1*sig_max]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("60Hz and 180Hz notch filter ECG");

subplot(4,1,4)
plot(time,signal2)
time = linspace(0,totaltime,m);
axis([0,totaltime,1.1*sig_min,1.1*sig_max]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("60Hz, 180Hz and 300Hz notch filter ECG");

%% 7)Resample
%The highest useful frequency:
f_useful = 50;  %If this is <50, the entire signal is ruined
%The corresponding new sample frequency:
fs_new = 2 * f_useful;

%Find decimation and interpolation values:
proportion = fs_new / fs;
[I,D] = findInterpolationDecimation(proportion,3);
I = double(I);
D = double(D);

%For any resampling, first do the interpolation
resampledSignal = upsample(signal,I);

%Next, use a low pass filter to prevent aliasing
%1) Calculated stopband frequencies for I and D:
fs_I = fs * I;          %Sample frequency after interpolation: Fy = Fx * I
fn_I = fs_I / (2*I);    %Stopband frequency of the new signal = Fy / (2I)
fs_D = fs_I / D;        %Sample frequency after decimation: Fy = Fx / D
fn_D = fs_D / 2;            %Stoband frequency of the new signal = Fy / 2

%The required minimal stopband frequency then is:
fstop = min(fn_I,fn_D);
df = 5;                %Transition band width in Hz
                        %This gives a cutoff frequency of:
                        %(this is the case for a transition band of 10Hz)
wc = (fstop + (fstop - df)) / fs_I * pi;
%wc = 0.25*pi;
%The low pass filters needs to be a IIR filter.
%For the formula, please refer to the lab report:
Ts_I = 1 / fs_I;
Ohm = (2 / Ts_I) * tan(wc / 2);
b = [Ohm * Ts_I, Ohm * Ts_I];
a = [Ohm * Ts_I + 2, Ohm * Ts_I - 2];
fvtool(b,a);
freqz(b,a);
resampledSignal = downsample(resampledSignal,D);

figure
subplot(2,1,1)
time = linspace(0,totaltime,m);
plot(time,signal)
axis([0,totaltime,1.1*sig_min,1.1*sig_max]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG");

subplot(2,1,2)
totalTime2 = fs_new * length(resampledSignal);
time2 = linspace(0,totalTime2,length(resampledSignal));
plot(time2,resampledSignal)
axis([0,totalTime2,1.1*min(resampledSignal),1.1*max(resampledSignal)]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG");

n2=2^nextpow2(length(resampledSignal)); %Efficiently calculate FFT using N = power of 2 
X_res = fft(resampledSignal,n2)/n2; %MATLAB requires this get the correct amplitude
X_res = abs(X_res); %Absolute value of the fft

f2 = fs_new*(0:(n2/2))/n2;
X_res = X_res(1:n2/2+1);
figure
plot(f2,X_res)