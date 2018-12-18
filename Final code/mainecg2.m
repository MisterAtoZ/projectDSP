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
axis([0,totaltime*0.05,1.1*sig_min,1.1*sig_max]);
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
%Rescaling based on reccomendations from
%https://nl.mathworks.com/help/matlab/ref/fft.html
f = fs*(0:(n/2))/n;
X_plot = X_abs(1:n/2+1);

%Plot the frequency spectrum
plot(f,X_plot)
title("Frequency spectrum of the ECG2");
legend("Absolute value");
xlabel("Frequentie (Hz)");


%% BP Filter:
%We need 2 cutoff frequencies: wc1 en wc2:
f_pass1 = 30;
f_stop1 = 35;
f_pass2 = 0;
f_stop2 = 4;
w1s=f_stop1/fn*pi;     %Stopband1: [w1s,pi]
w1p=f_pass1/fn*pi;     %Passband1: [0,w1p]
w2p=f_pass2/fn*pi;     %Passband2: [0,w2p]
w2s=f_stop2/fn*pi;     %Stopband2: [w2s,pi]
% The cut-off frequency is halfway the transition band
wc1=(w1p+w1s)/2;
wc2=(w2p+w2s)/2;
%The transition band the filter needs to be configured on,
%is the smallest of the two transition bands
dw=min(w1s-w1p,w2s-w2p);
%We choose a stopband attenuation of 60dB
As = 60;

%We use a Kaiser window
M=(As-7.95)/(2.285*dw) + 1;
M = roundToNextOddInteger(M);
if As >= 50
    beta = 0.1102*(As-8.7);
elseif (As > 21) && (As < 50)
    beta = 0.5842*(As-21)^0.4 + 0.07886*(As-21);
else
    error('Error: the chosen As=%d is smaller than 22\n',As);
end
W = kaiser(M,beta);

alfa = (M-1)/2;
b_lp1 = wc1 / pi * sinc(wc1 / pi * (-alfa:alfa));
b_lp2 = wc2 / pi * sinc(wc2 / pi * (-alfa:alfa));

b_bp = (b_lp1 - b_lp2).*W';
fvtool(b_bp);

figure
signal_bp = filter(b_bp,1,signal);
plot(time,signal_bp)
hold on
plot(time,signal)
axis([0,totaltime*0.05,1.1*min(min(signal),min(signal_bp)),1.1*max(max(signal),max(signal_bp))]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2 before and after bandpass filter");
legend("ECG2 after bandpass filter","ECG2 before bandpass filter");

figure
X_bp = abs(fft(signal_bp,n) / n);
X_bp = X_bp(1:n/2+1);
plot(f,X_plot)
hold on
plot(f,X_bp)
hold off
xlabel("Frequency in Hz",'fontSize',14)
ylabel("Normalized magnitude",'fontSize',14)
title("Frequency spectrum before and after bandpass")
legend("Before bandpass","After bandpass",'fontSize',14)
figure
subplot(2,1,1)
plot(time,signal)
axis([0,totaltime*0.05,1.1*min(signal),1.1*max(signal)]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2 before bandpass");

subplot(2,1,2)
plot(time,signal_bp)
axis([0,totaltime*0.05,1.1*min(signal_bp),1.1*max(signal_bp)]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2 after bandpass");

figure
subplot(2,1,1)
plot(time,signal)
axis([0,totaltime,1.1*min(signal),1.1*max(signal)]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2 before bandpass");

subplot(2,1,2)
plot(time,signal_bp)
axis([0,totaltime,1.1*min(signal_bp),1.1*max(signal_bp)]);
xlabel("Time in s");
ylabel("Signal amplitude");
title("ECG2 after bandpass");