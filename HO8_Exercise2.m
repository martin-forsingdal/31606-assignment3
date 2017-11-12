% Exercise 2
% Clear all variables and command window 
clear all;
close all;
clc;

%% Test the combination of the filters

% Set the filter order
order=1000;
% Set the normalized sampling frequency
fs = 2;
% Compute the length of the filter
N0 = order + 1;
% Set the cut-off frequencies
w_l = 0.011338;
w_b1 = 0.090703;
w_b2 = 0.181406;
w_b3 = 0.272109;
w_b4 = 0.725624;
% Set the offset frequency
offset = 0.004/5;
% Create the individual filters 
[h_l,~] = FIR_lowpass(w_l-offset,-3,0,order);
[h_b1,~] = FIR_bandpass([w_l+offset,w_b1-offset],-1.5,0,order);
[h_b2,~] = FIR_bandpass([w_b1+offset,w_b2-offset],1,0,order);
[h_b3,~] = FIR_bandpass([w_b2+offset,w_b3-offset],0.5,0,order);
[h_b4,~] = FIR_bandpass([w_b3+offset,w_b4-offset],-0.3,0,order);

% Compute the frequency response by Fourier transform
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

% Plot the filter on top of each other
figure(1);
subplot(3,1,[1,2]);
% Plot the combined filter
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H))),'k','LineWidth',6);
hold on;
grid on;
% Plot the lowpass filter
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_l))),'c','LineWidth',2);
% Plot the first bandpass filter
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_b1))),'r','LineWidth',2);
% Plot the second bandpass filter
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_b2))),'g','LineWidth',2);
% Plot the third bandpass filter
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_b3))),'m','LineWidth',2);
% Plot the fourth bandpass filter
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_b4))),'y','LineWidth',2);
xlabel('Normalized frequency [\times \pi rad/sample]','FontSize',12);
ylabel('Gain [dB]','FontSize',15);
legend({'Combined Bandpass','Lowpass','First bandpass',...
    'Second bandpass','Third bandpass','Fourth bandpass'}...
    ,'FontSize',12);
title('Magnitude','FontSize',15);
xlim([0 1])
hk = ifft(H);
hold off
subplot(3,1,3);
% Plot the phase of the combined filter
plot((-order/2:order/2)/order*fs,fftshift(rad2deg(unwrap(angle(H)))),'k',...
    'LineWidth',1);
xlabel('Normalized frequency [\times \pi rad/sample]','FontSize',12);
ylabel('Phase [deg]','FontSize',12);
grid on;
xlim([0 1]);
title('Phase','FontSize',15);
% Export the figure
% hgexport(gcf,'filter_responses');
%% Test the FIR_eq function
% Create a filter with the attenuations given by the inputs to FIR_eq
[wk, H] = FIR_eq(-3,-1.5,1,0.5,0.3);
% Plot the resulting filter
figure(2);
subplot(3,1,[1 2]);
plot(wk,fftshift(mag2db(abs(H))),'k','LineWidth',2);
xlabel('Normalized frequency [\times \pi rad/sample]','FontSize',12);
ylabel('Gain [dB]','FontSize',12);
grid on;
xlim([0 1]);
title('Magnitude','FontSize',15);
subplot(3,1,3);
plot(wk,fftshift(rad2deg(unwrap(angle(H)))),'k','LineWidth',2)
xlabel('Normalized frequency [\times \pi rad/sample]','FontSize',12);
ylabel('Phase [deg]','FontSize',12);
grid on;
xlim([0 1]);
title('Phase','FontSize',15);
% Export the figure
% hgexport(gcf,'equalizer_filter');