%% Exercise 1.2
% Clear all variables and command window
clear all;
close all;
clc;

%% Determine sampling frequency
% The known specifications are:
% Passband gain
Rp_dB=0;
Rp=db2mag(Rp_dB);
% Passband gain exists up to the passband frequency (in Hz)
fp = 5e3;
% Normalized cufoff frequency is
wn_norm=1/3;
% Number of samples in impulse response
N0=601;
% Determine sampling freqeuncy fs.
% wn_norm*fs/2=fp
% Isolate fs, and we get
fs=2*fp/wn_norm; 
%% Construct filter using frequency equivalence method
close all;
% Initialize the vector containing the sampled frequency response of
% ideal lowpass filter
H = ones(1,N0);
% Set the zeros
H(102:end-100) = 0;
% Frequency resolution is
F0=fs/N0;
% Define frequency vector.
freq=0:2*pi*F0:2*pi*(fs-F0);
% Plot the magnitude of the frequency response
figure('Name', 'Spectrum', 'Position', [800 200 700 500])
subplot(2,1,1)
stem(freq-pi*fs, fftshift(abs(H)), 'r.')
xlabel('\omega [rad/s]');
ylabel('|H| [linear scale]')
ylim([-0.25, 1.25])
grid on
set(gca, 'Fontsize', 14)
% Plot the phase of the frequency response
subplot(2,1,2)
stem(freq-pi*fs, fftshift(rad2deg(angle(H))), 'r.')
xlabel('\omega [rad/s]');
ylabel('\angle H [degrees]')
ylim([-190, 190])
grid on
set(gca, 'Fontsize', 14)
% Export
hgexport(gcf, 'handson8_ideal_lowpass.eps')
%%
% Define Hr, the version of the frequency response that is phase shifted
% by a linear phase shift: 
% exp(j*r*pi*(N0-1)/No, where r is {0, 1, ..., N0-1}.
% This version has a causal impulse response.
Hr = zeros(1,N0);
% Compute the elements using 12.97 in Lathi
for i = 1:N0
    % Set the r-variable
    r = i-1;
    % Compute the elements
    Hr(i) = H(i)*exp(j*r*pi*(N0-1)/N0);
end
% Define the vector containing all sample numbers from 0 to N0-1
k = 0:N0-1;
% Find the impulse response by inverse Fourier transform (IDFT)
h = real(ifft(Hr));
% Check the largest appearing imaginary part
max(abs(imag(ifft(Hr))));
% it is in the order of 10^-14, so these imaginary parts due to 
% numerical errors are neglected
%% Window impulse response
% window the impulse response using rect functions, so impulse responses
% with 101, 21 and 11 nonzero samples are obtained
h_101=h.*generate_rect(N0, 101);
h_21=h.*generate_rect(N0, 21);
h_11=h.*generate_rect(N0, 11);
%% Plot the impulse responses
close all;
figure('Name', 'Impulse response', 'Position', [800 200 800 800])
% plot original 601 non-zero element impulse response
subplot(4,1,1)
stem(k, h, 'k.');
xlabel('k [sample number]');
ylabel('Amplitude [a.u.]');
ylim([-.1, .4])
legend('601 non-zero elements')
set(gca, 'Fontsize', 14)
grid on
% Plot 101 non-zero element impulse response
subplot(4,1,2)
stem(k, h_101, 'b.');
xlabel('k [sample number]');
ylabel('Amplitude [a.u.]');
ylim([-.1, .4])
legend('101 non-zero elements')
set(gca, 'Fontsize', 14)
grid on
% Plot 21 non-zero element impulse response
subplot(4,1,3)
stem(k, h_21, 'r.');
xlabel('k [sample number]');
ylabel('Amplitude [a.u.]');
ylim([-.1, .4])
legend('21 non-zero elements')
set(gca, 'Fontsize', 14)
grid on
% Plot 11 non-zero element impulse response
subplot(4,1,4)
stem(k, h_11, 'm.');
xlabel('k [sample number]');
ylabel('Amplitude [a.u.]');
ylim([-.1, .4])
legend('11 non-zero elements')
set(gca, 'Fontsize', 14)
grid on
% Export
hgexport(gcf, 'handson8_windowing.eps')
%% Calculate spectrums
[H, freq]=spectrum_maker(h, fs*2*pi);
[H_101, ~]=spectrum_maker(h_101, fs*2*pi);
[H_21, ~]=spectrum_maker(h_21, fs*2*pi);
[H_11, ~]=spectrum_maker(h_11, fs*2*pi);
%% Plot the frequency responses on top of each other
close all;
figure('Name', 'frequency responses', 'Position', [700, 100, 900, 850])
% Plot the magnitude response
subplot(2,1,1)
P_m601=plot(freq-pi*fs, fftshift(mag2db(abs(H))), 'k-');
hold on
P_m101=plot(freq-pi*fs, fftshift(mag2db(abs(H_101))), 'b-');
P_m21=plot(freq-pi*fs, fftshift(mag2db(abs(H_21))), 'r-');
P_m11=plot(freq-pi*fs, fftshift(mag2db(abs(H_11))), 'm-');
xlabel('\omega [rad/s]')
ylabel('|H| [dB]')
ylim([-120, 0])
legend('DFT of h[k]', 'DFT of 101 non-zero element h[k]',...
    'DFT of 21 non-zero element h[k]',...
    'DFT of 11 non-zero element h[k]',...
    'Location', 'northeast') 
title('Magnitude response')
grid on
set(gca, 'Fontsize', 14)
lw=1.2;
set(P_m601, 'Linewidth', lw)
set(P_m101, 'Linewidth', lw)
set(P_m21, 'Linewidth', lw)
set(P_m11, 'Linewidth', lw)
% Plot the phase of the frequency response
subplot(2,1,2)
P_p601=plot(freq-pi*fs, fftshift(rad2deg(unwrap(angle(H)))), 'k-');
hold on
P_p101=plot(freq-pi*fs, fftshift(rad2deg(unwrap(angle(H_101)))), 'b-');
P_p21=plot(freq-pi*fs, fftshift(rad2deg(unwrap(angle(H_21)))), 'r-');
P_p11=plot(freq-pi*fs, fftshift(rad2deg(unwrap(angle(H_11)))), 'm-');
xlabel('\omega [rad/s]')
ylabel('\angle H [degrees]')
title('Phase response')
grid on
legend('DFT of h[k]', 'DFT of 101 non-zero element h[k]',...
    'DFT of 21 non-zero element h[k]',...
    'DFT of 11 non-zero element h[k]',...
    'Location', 'northeast') 
set(gca, 'Fontsize', 14)
set(P_p601, 'Linewidth', lw)
set(P_p101, 'Linewidth', lw)
set(P_p21, 'Linewidth', lw)
set(P_p11, 'Linewidth', lw)
% Export
hgexport(gcf, 'handson8_windowed_specs.eps')


