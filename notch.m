function [y] = notch(signal, fs, f0, notchWidth)
% tabel = load('ecg.mat');
% signal = tabel.ecg;

% fs = 20000;             %#sampling rate
% f0 = 50;                %#notch frequency
fn = fs/2;              %#Nyquist frequency
freqRatio = f0/fn;      %#ratio of notch freq. to Nyquist freq.

%Hoe breder deze wordt genomen, hoe meer frequenties errond ook worden
%verzwakt
%notchWidth = 0.1;       %#width of the notch

%Compute zeros
notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

%#Compute poles
notchPoles = (1-notchWidth) * notchZeros;

% figure;
% zplane(notchZeros.', notchPoles.');

b = poly( notchZeros ); %# Get moving average filter coefficients
a = poly( notchPoles ); %# Get autoregressive filter coefficients

% figure;
% freqz(b,a,32000,fs)

%#filter signal x
y = filter(b,a,signal);
