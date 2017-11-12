function [hk, Hr] = FIR_bandpass(Wc,Rp,Rs,order)
% Wc is vector of 2 normalized frequencies, nomalized with Nyquist
% frequency
% Rs is passband attenuation
% Rp is stopband attenuation
% order is the order of the filter, the length of the impulse response will
% be order + 1

N0=order+1;

lowIndex = ceil(N0/2*Wc(1));
highIndex = ceil(N0/2*Wc(2));

H = ones(1,N0) * db2mag(-Rs);
H(lowIndex:highIndex) = 1 * db2mag(-Rp);
H(N0-highIndex+1:N0-lowIndex+1) = 1 * db2mag(-Rp);

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