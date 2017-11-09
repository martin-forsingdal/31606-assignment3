%% Exercise 1.1
% Clear all variables and command window
clear all;
close all;
clc;

%% Construct Butterworth filter
Wp = [0.2 0.3];
Ws = [0.1 0.4];
Rp = 2;
Rs = 100;
[n, Wn] = buttord(Wp,Ws,Rp,Rs);
[B, A] = butter(n,Wn,'bandpass');
[h, w] = freqz(B,A,1e3);

%% Plot the frequency response using freqz
figure(1)
freqz(B,A);

%% Plot the frequency response using output from freqz
close all;
figure(1)
plot(w/pi,mag2db(abs(h)),'r-','LineWidth',1.5);
xlabel('Normalized frequency [\times \pi rad/sample]','FontSize',15);
ylabel('Amplitude [dB]','FontSize',15);
grid on;
titl = sprintf(['Frequency response of Butterworth bandpass filter'...
    ' of order %d'],2*n);
title(titl,'FontSize',20);

%% Compute the impulse response
[u t] = impz(B,A,1e5);

% Plot impulse response
maxu = max(abs(u));
index = find(abs(u)>0.1*maxu,1,'last');
figure(2);
plot(t(1:index+1),u(1:index+1),'b-','LineWidth',1.5);
xlabel('Time [s]','FontSize',15);
ylabel('Impulse value [a.u.]','FontSize',15);
grid on;
title('Plot of impulse response cut to its effective length','FontSize',20);

%% Retun to frequency domain
hand = freqz(u);
