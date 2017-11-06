%hands_on_08
clc
clear
%% 1.1
close all;
% Create filter coefficients for Butterworth bandpass filter 
% with the given specifications
wp=[.2 .3];
ws=[.1 .4];
Rp=2;
Rs=100;
[N, wn]=buttord(wp, ws, Rp, Rs);
[B, A]=butter(N, wn);
% Compute the frequency response
[H, w]=freqz(B, A);
H=H';
w=w'./pi;
% Plot the magnitude and the pahse of the frequency response
figure('Name', 'Frequncy response', 'position', [800 300 700 450]);
subplot(2,1,1)
plot(w, mag2db(abs(H)), 'b-')
ylim([-200, 20])
xlim([0 1])
set(gca, 'Fontsize', 14)
% Add vertical lines to illustrate filter specifiations
y_lims=get(gca, 'YLim');
line([wp(1), wp(1)], [Rp/2+ 1/3*y_lims(1), Rp/2],...
    'color', [.5, .5, .5])
line([wp(2), wp(2)], [Rp/2+ 1/3*y_lims(1), Rp/2],...
    'color', [.5, .5, .5])
line([ws(1), ws(1)], [-Rs, Rp/2-25],...
    'color', [.5, .5, .5])
line([ws(2), ws(2)], [-Rs, Rp/2-25],...
    'color', [.5, .5, .5])
% Add horisontal lines to illustrate filter specifications
x_lims=get(gca, 'Xlim');
line([wp(1), wp(2)], [Rp/2, Rp/2],...
    'color', [.5, .5, .5])
line([x_lims(1), ws(1)], [-Rs, -Rs],...
    'color', [.5, .5, .5])
line([ws(2), ws(2)+1/10], [-Rs, -Rs],...
    'color', [.5, .5, .5])
% The order of the filter is the double of the value returned
% by buttord (see documentation for butter function)
title(string(2*N) + 'th order Butterworth filter')
xlabel('Normalized frequency')
ylabel('|H| [dB]')
subplot(2,1,2)
plot(w, rad2deg(angle(H)), 'r-')
xlabel('Normalized frequency')
ylabel('Angle [degrees]')
set(gca, 'Fontsize', 14)
grid on
% Compute the impulse response
[h, t]=impz(B, A);
% Plot the impulse response
figure('Name', 'Impulse response', 'position', [800 300 500 330]);
stem(t, h, 'r.')
xlabel('k')
ylabel('A a.u.')
title(string(2*N) + 'th order Butterworth filter')
grid on
set(gca, 'Fontsize', 13)
grid on
% pole/zero location 
figure('Name', 'Pole and zero location', 'position', [800 200 400 400]);
zplane(B, A)
set(gca, 'Fontsize', 14)
grid on
%% Cut impulse response at effective length
close all;
% Truncate impulse response when amplitude is no longer higher than
% 10% of the maximum value, to obtain the effective impulse response. 
% Find the last index where this happens
b=abs(h)>0.1*max(abs(h));
i=find(b, 1, 'last');
% Define effective response as the original response up to the found
% index i
h_eff=h(1:i);
%% Window the effective impulse response
% Define new vectors, which are the effective impulse response cutted at
% 75%, 60 %, 50% and 10 % of the length of the effective impulse response
h_eff75=h_eff(1:floor(.75*(length(h_eff)-1)));
h_eff60=h_eff(1:floor(.6*(length(h_eff)-1)));
h_eff40=h_eff(1:floor(.4*(length(h_eff)-1)));
h_eff10=h_eff(1:floor(.1*(length(h_eff)-1)));
% Transform to the Fourier domain to obtain the filters of all versions of
% the effective impulse response. Use sampling frequency fs=2 to obtain 
% normalized spectrum
fs=2;
[H_eff, freq_eff]=spectrum_maker(h_eff, fs);
[H_eff75, freq_eff75]=spectrum_maker(h_eff75, fs);
[H_eff60, freq_eff60]=spectrum_maker(h_eff60, fs);
[H_eff40, freq_eff40]=spectrum_maker(h_eff40, fs);
[H_eff10, freq_eff10]=spectrum_maker(h_eff10, fs);

% Plot the original effective impulse response and all cutted versions
figure('Name', 'Cutted effective impulse responses',...
    'Position', [800 200 900 600])
% The magnitude of the spectrums
subplot(2,1,1)
pmag=plot(freq_eff-.5*fs, mag2db(abs(H_eff)), '-', 'color', [.5, .5, .5]);
hold on
pmag75=plot(freq_eff75-.5*fs, mag2db(abs(H_eff75)), 'g-');
pmag60=plot(freq_eff60-.5*fs, mag2db(abs(H_eff60)), 'b-');
pmag40=plot(freq_eff40-.5*fs, mag2db(abs(H_eff40)), 'r-');
pmag10=plot(freq_eff10-.5*fs, mag2db(abs(H_eff10)), 'm-');
xlabel('Normalized frequency');
ylabel('|H| [dB]')
xlim([0, 1]);
ylim([-120, -35])
legend('|H_{eff}|', '75% of |H_{eff}|', '60% of |H_{eff}|',...
    '40% of |H_{eff}|', '10% of |H_{eff}|',...
    'Location', 'NorthEastOutside')
set(gca, 'Fontsize', 14)
set(pmag, 'Linewidth', 2)
set(pmag75, 'Linewidth', 2)
set(pmag60, 'Linewidth', 2)
set(pmag40, 'Linewidth', 2)
set(pmag10, 'Linewidth', 2)
grid on
% The phase of the spectrums
subplot(2,1,2)
pmag=plot(freq_eff-.5*fs, rad2deg(angle(H_eff)), '-', 'color', [.5, .5, .5]);
hold on
pmag75=plot(freq_eff75-.5*fs, rad2deg(angle(H_eff75)), 'g-');
pmag60=plot(freq_eff60-.5*fs, rad2deg(angle(H_eff60)), 'b-');
pmag40=plot(freq_eff40-.5*fs, rad2deg(angle(H_eff40)), 'r-');
pmag10=plot(freq_eff10-.5*fs, rad2deg(angle(H_eff10)), 'm-');
xlabel('Normalized frequency');
ylabel('Angle [degrees]')
legend('\angle H_{eff}', '75% of \angle H_{eff}', '60% of \angle H_{eff}',...
    '40% of \angle H_{eff}', '10% of \angle H_{eff}',...
    'Location', 'NorthEastOutside')
set(gca, 'Fontsize', 14)
set(pmag, 'Linewidth', 1)
set(pmag75, 'Linewidth', 1)
set(pmag60, 'Linewidth', 1)
set(pmag40, 'Linewidth', 1)
set(pmag10, 'Linewidth', 1)
grid on
%% 1.2 shorter and shorter
clc
clear
%%
close all;
% The known specifications are:
% Passband gain
Rp_dB=0;
Rp=db2mag(Rp_dB);
% Passband gain exists up to the passband frequency (in rad/s)
wp=2*pi*5*10^3;
% Normalized cufoff frequency is
wn_norm=1/3;
% Number of samples in impulse response
N_0=601;
% Calculate fs
fs=2*pi*N_0*wp/(1/sqrt(2)-Rp+2*pi^2*N_0*wn_norm*Rp);
% The frequency resolution and the frequency vector is then
F_0=fs/N_0;
freq=0:2*pi*F_0:2*pi*fs-2*pi*F_0;
% Create sampled lowpass filter with linear phase
% DC component
H=ones(1, N_0)*100;
% Define the positive frequencies. Find elements less than or 
% equal to passband frequency
i_pos=1:ceil(length(freq)/2);
b=freq(i_pos)<=wp;
% Set the first occuring zero to one, since Rp must excist
% throughout the WHOLE passband
i_firstzero=find(b==0, 1);
b(i_firstzero)=1;
% Define the complex spectrum for positive frequencies
H(1:(N_0-1)/2+1)=b*Rp;
H((N_0-1)/2+2:N_0)=fliplr(conj(b(2:length(b))*Rp));
% Plot the spectrum of the sampled ideal lowpass filter
figure('Name', 'Spectrum', 'Position', [800 200 700 500])
subplot(2,1,1)
stem(freq-pi*fs, fftshift(abs(H)), 'r.')
xlabel('\omega [rad/s]');
ylabel('|H| [dB]')
ylim([0, 2])
grid on
set(gca, 'Fontsize', 14)
subplot(2,1,2)
stem(freq-pi*fs, fftshift(rad2deg(angle(H))), 'r.')
xlabel('\omega [rad/s]');
ylabel('\angle H [degrees]')
ylim([-190, 190])
grid on
set(gca, 'Fontsize', 14)
% Transform to time domain to obtain impulse response
h_c=ifft(H);
h=real(h_c);
h_im=imag(h_c);
max(abs(h_im));
% Imaginary part of impulse response is zero as desired.
% Create time vector
t=0:1/fs:N_0/fs-1/fs;
% Plot the impulse response
figure('Name', 'Impulse response', 'Position', [800 200 800 300])
stem(t, h, 'b.')
grid on
xlabel('t [s]')
ylabel('|A| a.u.')
set(gca, 'Fontsize', 14)