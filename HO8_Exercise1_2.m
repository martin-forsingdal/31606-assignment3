%% Exercise 1.2
% Clear all variables and command window
clear all;
close all;
clc;

%% Determine sampling frequency
fp = 5e3;
N0 = 601;
fs = 6*fp;

%% Construct filter using frequency equivalence method
% Initialize transferfunction
H = ones(1,N0);
% Set the zeros
H(102:end-100) = 0;
% Plot the magnitude of the transfer function
figure(1);
stem(-300:300,fftshift(H),'r');
xlabel('Sample');
ylabel('|H_r|');
title('Plot of the magnitude of the desired transfer function');
grid on;
% Initialise vector for the frequency response
Hr = zeros(1,N0);
% Compute the elements using 12.97 in Lathi
for i = 1:N0
    % Set the r-variable
    r = i-1;
    % Compute the elements
    Hr(i) = H(i)*exp(j*r*pi*(N0-1)/N0);
end
% Set number of samples
k = 0:N0-1;
% Find the impulseresponse by inverse Fourier transform
hk = real(ifft(Hr));

%% Plot the impulse response
figure(2);
stem(k,hk);
xlabel('k');
ylabel('h[k]');
title('Plot of the impulse response of the filter');
% Set the number of points for which the filter transfer function should be
% computed
M = 10000;
% Define longer impulse response
hE = [hk zeros(1,M-N0)];
% Compute the transfer function
HE = fft(hE);

%% Plot the transfer function
figure(3)
subplot(2,1,1);
% Compute the frequency values
r = 0:M-1;
W = r.*2*pi/M;
% Plot the magnitude
plot((W/(pi)-1)*fs/2,fftshift(abs((HE))));
xlabel('f in [Hz]');
ylabel('|H[e^{j\omegaT}]|');
grid on;
title('Plot of the magnitude of the transfer function');
% Plot the phase
subplot(2,1,2);
plot((W/(pi)-1)*fs/2,fftshift(180/pi*unwrap(angle((HE)))));
xlabel('f in [Hz]');
ylabel('\angle H[e^{j\omegaT}] in [deg]');
grid on;
title('Plot of the phase of the transfer function');
