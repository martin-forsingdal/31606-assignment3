function [hk, Hr] = FIR_lowpass(Wc,Rp,Rs,order)
% This function creates a FIR lowpass filter with linear phase based on
% the given cut-off frequency pass- and stopband gains, and order.
% Input:    Wc = Normalized cut-off frequency, nomalized to Nyquist
%           frequency
%           Rs = Passband attenuation
%           Rp = Stopband attenuation
%           order = The order of the filter, the length of the impulse 
%           response will be order + 1      
% Output:   hk = Impulse response of the filter
%           Hr = The filter coefficients

% Define the length of the filter
N0 = order + 1;

% Compute the index of the cut-off frequency
lowIndex = ceil(N0/2*Wc);

% Preallocate the magnitude vector and set the stopband attenuation
H = ones(1,N0)*db2mag(-Rs);
% Set the passband gain
H(1:lowIndex) = 1*db2mag(-Rp);
% Make the spectrum symmetric to get the linear phase
H(end-lowIndex+1:end) = 1*db2mag(-Rp);

% Preallocate the vector containing the filter coefficients
Hr = zeros(1,N0);
% Compute the elements using 12.97 in Lathi
for i = 1:N0
    % Set the r-variable
    r = i-1;
    % Compute the elements
    Hr(i) = H(i)*exp(j*r*pi*(N0-1)/N0);
end
% Find the impulseresponse by inverse Fourier transform
hk = real(ifft(Hr));
% Window the impulse response using a Hann window of length N0
hk = hk.*hann(N0)';

end
