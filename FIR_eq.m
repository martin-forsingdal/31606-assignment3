function [freq, H] = FIR_eq(R_l,R_b1,R_b2,R_b3,R_b4)

order=1000;
N0 = order + 1;
w_l = 0.011338;
w_b1 = 0.090703;
w_b2 = 0.181406;
w_b3 = 0.272109;
w_b4 = 0.725624;
offset = 0.004/5;
[h_l,~] = FIR_lowpass(w_l-offset,R_l,0,order);
[h_b1,~] = FIR_bandpass([w_l+offset,w_b1-offset],R_b1,0,order);
[h_b2,~] = FIR_bandpass([w_b1+offset,w_b2-offset],R_b2,0,order);
[h_b3,~] = FIR_bandpass([w_b2+offset,w_b3-offset],R_b3,0,order);
[h_b4,~] = FIR_bandpass([w_b3+offset,w_b4-offset],R_b4,0,order);

H_l = fft(h_l);
H_b1 = fft(h_b1);
H_b2 = fft(h_b2);
H_b3 = fft(h_b3);
H_b4 = fft(h_b4);

H = H_l.*H_b1;
H = H.*H_b2;
H = H.*H_b3;
H = H.*H_b4;

freq = -1:1/(order/2):1;

end