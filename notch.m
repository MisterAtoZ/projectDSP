function [y] = notch(signal, fs, f0, notchWidth)
fn = fs/2; %de maximum frequentie die kan worden gehaald worden zonder aliasing
freqRation = f0/fn;

%de juiste hoek zoeken van de nulpunten en polen
hoek = (pi/fn) * f0;
reelDeel = cos(hoek);
imagDeel = sin(hoek);

nulp1 = reelDeel + imagDeel*i
nulp2 = reelDeel - imagDeel*i
pool1 = notchWidth*nulp1;
pool2 = notchWidth*nulp2;

zplane(nulp1,pool1)
figure
zplane(nulp2,pool2)

%de TF van de filter wordt dan
%H = ((z-nulp1)*(z-nulp2))/((z-pool1)*(z-pool2));

%Deze filter toepassen op signal



% % tabel = load('ecg.mat');
% % signal = tabel.ecg;
% 
% % fs = 20000;             %#sampling rate
% % f0 = 50;                %#notch frequency
% fn = fs/2;              %#Nyquist frequency
% freqRatio = f0/fn;      %#ratio of notch freq. to Nyquist freq.
% 
% %Hoe breder deze wordt genomen, hoe meer frequenties errond ook worden
% %verzwakt
% %notchWidth = 0.1;       %#width of the notch
% 
% %Compute zeros
% notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];
% 
% %#Compute poles
% notchPoles = (1-notchWidth) * notchZeros;
% 
% % figure;
% % zplane(notchZeros.', notchPoles.');
% 
% b = poly( notchZeros ); %# Get moving average filter coefficients
% a = poly( notchPoles ); %# Get autoregressive filter coefficients
% 
% % figure;
% % freqz(b,a,32000,fs)
% 
% %#filter signal x
% y = filter(b,a,signal);
