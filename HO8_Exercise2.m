% Exercise 2
% Clear all variables and command window 
clear all;
close all;
clc;

%% Test
N0 = 302;
[h, Hr] = FIR_bandpass([0.2,0.6],-20,N0);
freqz(h.*hamming(N0)');
