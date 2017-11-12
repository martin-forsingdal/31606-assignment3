function [freq, H] = FIR_eq(R_l,R_b1,R_b2,R_b3,R_b4)
% This function creates an FIR equalize filter which is composed by a 
% lowpass filter and four bandpass filter. The FIR filter is of order 1000.
% The band ranges:
%           w_l = DC to 0.011338
%           w_b1 = 0.011338 to 0.090703
%           w_b2 = 0.090703 to 0.181406
%           w_b3 = 0.181406 to 0.272109 
%           w_b4 = 0.272109 to 0.725624
% The ranges above are normalized to a sampling frequency of 44.1kHz
% Input:    R_l = Attenuation of the lowpass range
%           R_b1 = Attenuation of the first band
%           R_b2 = Attenuation of the second band
%           R_b3 = Attenuation of the third band
%           R_b4 = Attenuation of the fourth band
% Ouput:    freq = Vector containing the normalized freuqencies. 
%           H = Vector containing the filter coefficients

% Define filter order
order=1000;
% Compute the length of the filter
N0 = order + 1;
% Set the normalized cut-off frequencies
w_l = 0.011338;
w_b1 = 0.090703;
w_b2 = 0.181406;
w_b3 = 0.272109;
w_b4 = 0.725624;
% Set the off-set frequency
offset = 0.004/5;
% Create the filters to be multiplied
[h_l,~] = FIR_lowpass(w_l-offset,R_l,0,order);
[h_b1,~] = FIR_bandpass([w_l+offset,w_b1-offset],R_b1,0,order);
[h_b2,~] = FIR_bandpass([w_b1+offset,w_b2-offset],R_b2,0,order);
[h_b3,~] = FIR_bandpass([w_b2+offset,w_b3-offset],R_b3,0,order);
[h_b4,~] = FIR_bandpass([w_b3+offset,w_b4-offset],R_b4,0,order);

% Get the frequency response by Fourier transform
H_l = fft(h_l);
H_b1 = fft(h_b1);
H_b2 = fft(h_b2);
H_b3 = fft(h_b3);
H_b4 = fft(h_b4);

% Multiply the filter responses to get the combined filter
H = H_l.*H_b1;
H = H.*H_b2;
H = H.*H_b3;
H = H.*H_b4;

% Define the nomalized frequency axis
freq = -1:1/(order/2):1;

end