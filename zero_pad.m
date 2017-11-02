function [time_vector_modified, signal_modified] = zero_pad (t, s, T)
% zero_pad fills a signal vector with zeros after the last element 
% of the signal
% function call
% [time_vector_modified, signal_modified] = zero_pad (t, s, T)
% t is the original time vector
% s is the original signal
% T is the desired duration of the zer padded signal

% Determine time between two consequent samples
T_0=t(length(t))-t(length(t)-1);
% Determine sampling frequency
f_s=1/T_0;
% Create time vector
time_vector_modified=[0:1/f_s:T-1/f_s];
N_tot=floor(f_s*T);
N_zeros_to_add=length(time_vector_modified)-length(s);
% Add the calculated number of zeros
signal_modified=[s, zeros(1,N_zeros_to_add)];


