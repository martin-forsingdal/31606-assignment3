%% Exercise 1.1
% Clear all variables and command window
clear all;
close all;
clc;

%% Construct Butterworth filter
% Create filter coefficients for Butterworth bandpass filter 
% with the given specifications
wp = [0.2 0.3];
ws = [0.1 0.4];
Rp = 2;
Rs = 100;
[n, wn] = buttord(wp,ws,Rp,Rs);
[B, A] = butter(n,wn,'bandpass');
% Compute the 1000 element frequency response
[H, w] = freqz(B,A,1e3);
H=H';
w=w'./pi;

%% Plot the frequency response
close all;
% Plot the magnitude and the pahse of the frequency response
figure('Name', 'Frequncy response', 'position', [800 300 700 450]);
subplot(2,1,1)
p_mag=plot(w, mag2db(abs(H)), 'b-');
ylim([-200, 20])
xlim([0 1])
set(gca, 'Fontsize', 14)
set(p_mag, 'Linewidth', 2)
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
title(string(2*n) + 'th order Butterworth filter')
xlabel('Normalized frequency')
ylabel('|H| [dB]')
legend('|H|', 'Specifications')
subplot(2,1,2)
plot(w, rad2deg(angle(H)), 'r-')
xlabel('Normalized frequency')
ylabel('\angle H [degrees]')
set(gca, 'Fontsize', 14);
grid on
legend('\angle H')
% Export
hgexport(gcf, 'handson8_spec_orig.eps')
%% pole/zero location 
figure('Name', 'Pole and zero location', 'position', [800 200 500 500]);
p_pz=zplane(B, A);
set(gca, 'Fontsize', 14)
set(p_pz, 'Linewidth', 1.2)
grid on
% Export
hgexport(gcf, 'handson8_pz.eps')
%% Impulse response
% Compute the 10000 element impulse response
[h, t] = impz(B,A,1e5);
% Cut impulse response at effective length.
% Truncate impulse response when amplitude is no longer higher than
% 10% of the maximum value, to obtain the effective impulse response. 
% Find the last index where this happens
maxu = max(abs(h));
i_last_element = find(abs(h)>0.1*maxu,1,'last');
% Define effective response as the original response up to the found
% index i
h_eff=h(1:i_last_element);
% Window the effective impulse response
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
%% Plot impulse responses and frequency responses
close all;
% plot impulse responses on top of each other
figure('Name', 'impulse responses', 'Position', [700, 100, 800, 300])
p_eff=stem(0:length(h_eff)-1, h_eff, 'k.');
hold on
p_eff75=stem(0:length(h_eff75)-1, h_eff75, 'r.');
p_eff60=stem(0:length(h_eff60)-1, h_eff60, 'm.');
p_eff40=stem(0:length(h_eff40)-1, h_eff40, 'b.');
p_eff10=stem(0:length(h_eff10)-1, h_eff10, 'g.');
xlabel('k [sample number]')
ylabel('|A| a.u.')
legend('h_{eff}[k]', '75% of h_{eff}[k]', '60% of h_{eff}[k]',...
    '40% of h_{eff}[k]','10% of h_{eff}[k]',...
    'Location', 'NorthEastOutside') 
set(gca, 'Fontsize', 14)
grid on
% Export
hgexport(gcf, 'handson8_impulse_responses.eps')
% plot frequency responses on top of each other.
figure('Name', 'frequency responses', 'Position', [700, 100, 900, 850])
% Plot the magnitude response
subplot(2,1,1)
P_eff=plot(freq_eff-1, fftshift(mag2db(abs(H_eff))), 'k-');
hold on
P_eff75=plot(freq_eff75-1, fftshift(mag2db(abs(H_eff75))), 'r-');
P_eff60=plot(freq_eff60-1, fftshift(mag2db(abs(H_eff60))), 'm-');
P_eff40=plot(freq_eff40-1, fftshift(mag2db(abs(H_eff40))), 'b-');
P_eff10=plot(freq_eff10-1, fftshift(mag2db(abs(H_eff10))), 'g-');
xlabel('Normalized frequency')
ylabel('|H| [dB]')
title('Magnitude response')
grid on
legend('DFT of h_{eff}[k]', 'DFT of 75% of h_{eff}[k]',...
    'DFT of 60% of h_{eff}[k]',...
    'DFT of 40% of h_{eff}[k]','DFT of 10% of h_{eff}[k]',...
    'Location', 'NorthEastOutside') 
set(gca, 'Fontsize', 14)
lw=1.2;
set(P_eff, 'Linewidth', lw)
set(P_eff75, 'Linewidth', lw)
set(P_eff60, 'Linewidth', lw)
set(P_eff40, 'Linewidth', lw)
set(P_eff10, 'Linewidth', lw)
% Plot the phase of the frequency response
subplot(2,1,2)
P_eff=plot(freq_eff-1, fftshift(rad2deg(angle(H_eff))), 'k-');
hold on
P_eff75=plot(freq_eff75-1, fftshift(rad2deg(angle(H_eff75))), 'r-');
P_eff60=plot(freq_eff60-1, fftshift(rad2deg(angle(H_eff60))), 'm-');
P_eff40=plot(freq_eff40-1, fftshift(rad2deg(angle(H_eff40))), 'b-');
P_eff10=plot(freq_eff10-1, fftshift(rad2deg(angle(H_eff10))), 'g-');
xlabel('Normalized frequency')
ylabel('\angle H [degrees]')
title('Phase response')
grid on
legend('DFT of h_{eff}[k]', 'DFT of 75% of h_{eff}[k]',...
    'DFT of 60% of h_{eff}[k]',...
    'DFT of 40% of h_{eff}[k]','DFT of 10% of h_{eff}[k]',...
    'Location', 'NorthEastOutside') 
set(gca, 'Fontsize', 14)
lw=1.2;
set(P_eff, 'Linewidth', 1)
set(P_eff75, 'Linewidth', 1)
set(P_eff60, 'Linewidth', 1)
set(P_eff40, 'Linewidth', 1)
set(P_eff10, 'Linewidth', 1)
% Export
hgexport(gcf, 'handson8_spec_cutted.eps')