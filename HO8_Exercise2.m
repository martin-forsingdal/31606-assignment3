% Exercise 2
% Clear all variables and command window 
clear all;
close all;
clc;

%% Test
order=300;
[h, Hr] = FIR_bandpass([0.2,0.6],-10,20,order);
freqz(h.*hann(N0)');