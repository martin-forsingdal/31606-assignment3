% Exercise 2
% Clear all variables and command window 
clear all;
close all;
clc;

%% Test

[wk, H] = FIR_eq(-3,-1.5,1,0.5,0.3);

figure(1);
plot(wk,fftshift(mag2db(abs(H))),'k','LineWidth',2);
xlim([0 1]);

%% Test FIR_lowpass

order=1000;
fs = 2;
N0 = order + 1;
w_l = 0.011338;
w_b1 = 0.090703;
w_b2 = 0.181406;
w_b3 = 0.272109;
w_b4 = 0.725624;
offset = 0.004/5;
[h_l,~] = FIR_lowpass(w_l-offset,-3,0,order);
[h_b1,~] = FIR_bandpass([w_l+offset,w_b1-offset],-1.5,0,order);
[h_b2,~] = FIR_bandpass([w_b1+offset,w_b2-offset],1,0,order);
[h_b3,~] = FIR_bandpass([w_b2+offset,w_b3-offset],0.5,0,order);
[h_b4,~] = FIR_bandpass([w_b3+offset,w_b4-offset],-0.3,0,order);

H_l = fft(h_l);
H_b1 = fft(h_b1);
H_b2 = fft(h_b2);
H_b3 = fft(h_b3);
H_b4 = fft(h_b4);

H = H_l.*H_b1;
H = H.*H_b2;
H = H.*H_b3;
H = H.*H_b4;



figure(1);
subplot(3,1,[1,2]);
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H))),'k','LineWidth',6);
hold on;
grid on;
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_l))),'c','LineWidth',2);
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_b1))),'r','LineWidth',2);
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_b2))),'g','LineWidth',2);
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_b3))),'m','LineWidth',2);
plot((-order/2:order/2)/order*fs,fftshift(mag2db(abs(H_b4))),'y','LineWidth',2);
xlabel('Normalized frequency [\times \pi rad/sample]','FontSize',15);
ylabel('Gain [dB]','FontSize',15);
legend({'Combined Bandpass','Lowpass [0,0.1]','Bandpass [0.1,0.3]',...
    'Bandpass [0.3,0.497]','Bandpass [0.502,0.69]','Bandpass [0.703,0.9]'}...
    ,'FontSize',15);
xlim([0 1])
hk = ifft(H);
hold off
subplot(3,1,3);
plot((-order/2:order/2)/order*fs,fftshift(rad2deg(unwrap(angle(H)))),'k',...
    'LineWidth',1);
xlabel('Normalized frequency [\times \pi rad/sample]','FontSize',15);
ylabel('Phase [deg]','FontSize',15);
grid on;
xlim([0 1]);

