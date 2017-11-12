%% init
clc
clear
close all
%% Read sound file
close all;
%filename = 'Steely Dan - Jack Of Speed.mp3';
filename = 'Aqua - Doctor Jones.mp3';
[soundsignal, fs] = audioread(filename);
testsignal=soundsignal(floor(1/18*end):floor(1/10*end));
dur = length(testsignal)/fs;
t = 0:1/fs:dur-1/fs;
%% Play the sound file
soundsc(testsignal,fs);

%% create eq filter
[wk, H]=FIR_eq(-10,-5,0,5,10);
testsignal_filter=real(filter(ifft(H),1,testsignal));
figure(1);
freqz(ifft(H));
%% play filtered signal
soundsc(testsignal_filter,fs);

%% plot frequency response of original and filtered signal
figure(2);
freqz(fft(testsignal)/length(testsignal));
hold on;
figure(3);
freqz(fft(testsignal_filter)/length(testsignal_filter));
hold off;