close all
tabel = load('ecg.mat');
signal = tabel.ecg;

fs = 1000;             %#sampling rate
f0 = 60;                %#notch frequency
fn = fs/2;              %de maximum frequentie die kan worden gehaald worden zonder aliasing
a = 0.9;
freqRatio = f0/fn;      %#ratio of notch freq. to Nyquist freq.

hoek = (pi/fn) * f0;

zeros = [1 -2*cos(hoek) 1];
poles = [1 -2*a*cos(hoek) a*a];

zplane(zeros,poles);

figure
m = length(signal);
ymax = max(signal);
ymin = min(signal);
Ts = 1/fs;
totaltime = Ts*m;
time = linspace(0,totaltime,m);
plot(time,signal)
axis([0,totaltime,1.1*ymin,1.1*ymax]);
xlabel("Time in s");

hold on

filteredSignal = filter(zeros,poles,signal);
plot(time,filteredSignal)

figure

subplot(3,1,1)
plot(time,signal)
axis([0,totaltime,1.1*ymin,1.1*ymax]);
xlabel("Time in s");

subplot(3,1,2)
plot(time,filteredSignal)
axis([0,totaltime,1.1*ymin,1.1*ymax]);
xlabel("Time in s");

%filter some more:
filteredSignal2 = filter(zeros,poles,filteredSignal);

subplot(3,1,3)
plot(time,filteredSignal2)
axis([0,totaltime,1.1*ymin,1.1*ymax]);
xlabel("Time in s");

%Frequency response of the filter:
fvtool(zeros,poles)
figure
%FFT of the signal after 
%try different angles:
X_old = fft(signal);
X_new = fft(filteredSignal);
subplot(4,1,1)
plot(abs(X_old))
%Belangrijk: We mogen maar kijken tot fs, anders krijg je enkel spiegeling
axis([0,fs,0,1.1*max(abs(X_old))]);
subplot(4,1,2)
plot(abs(X_new))
axis([0,fs,0,1.1*max(abs(X_new))]);

%Poging om nu de harmonische weg te filteren
subplot(4,1,3)
f0 = 3*60;
hoek = (pi/fn) * f0;
zeros = [1 -2*cos(hoek) 1];
poles = [1 -2*a*cos(hoek) a*a];

filteredSignal = filter(zeros,poles,filteredSignal);
X_new2 = fft(filteredSignal);
plot(abs(X_new2))
axis([0,fs*2.4,0,1.1*max(abs(X_new2))]);

subplot(4,1,4)
f0 = 5*60;
hoek = (pi/fn) * f0;
zeros = [1 -2*cos(hoek) 1];
poles = [1 -2*a*cos(hoek) a*a];

filteredSignal = filter(zeros,poles,filteredSignal);
X_new3 = fft(filteredSignal);
plot(abs(X_new3))
axis([0,fs,0,1.1*max(abs(X_new3))]);

%Kijk naar het tijdsdomein
figure
plot(time,filteredSignal)
hold on
plot(time, signal)
axis([0,totaltime,1.1*min(filteredSignal),1.1*max(filteredSignal)]);
xlabel("Time in s");

