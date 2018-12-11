%function [y] = notch(signal, fs, f0)
tabel = load('ecg.mat');
signal = tabel.ecg;

fs = 1000;             %#sampling rate
f0 = 60;                %#notch frequency
fn = fs/2;              %de maximum frequentie die kan worden gehaald worden zonder aliasing
notchWidth = 0.9;
freqRatio = f0/fn;      %#ratio of notch freq. to Nyquist freq.

%de juiste hoek zoeken van de nulpunten en polen
hoek = (pi/fn) * f0;
reelDeel = cos(hoek);
imagDeel = sin(hoek);

nulp1 = reelDeel + imagDeel*i; %hier is waarsch wel een simpelere manier voor
nulp2 = reelDeel - imagDeel*i;
pool1 = notchWidth*nulp1;
pool2 = notchWidth*nulp2;

zplane(nulp1,pool1)
figure
zplane(nulp2,pool2)

%de TF van de filter wordt dan
%H = ((z-nulp1)*(z-nulp2))/((z-pool1)*(z-pool2));

%Deze filter toepassen op signal


