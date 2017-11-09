% Exercise 2
% Clear all variables and command window 
clear all;
close all;
clc;

%% Test
order=100;
[h, Hr] = FIR_bandpass([0.2,0.6],-10,0,order);
figure(1);
freqz(h);
figure(2);
freqz(h.*hann(order+1)');
