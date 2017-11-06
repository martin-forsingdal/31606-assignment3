function [Y, freq] = spectrum_maker(signal, fs)
% Function computes complex-valued Fourier-spectrum of signal "signal" with
% sampling frequency fs.
% Function call:
% [Y, freq] = spectrum_maker(signal, fs)
% Input: signal to transform, "signal", sampling frequency of signal, fs:
% (signal, fs).
% Output: Complex-valued spectrum, Y, frequency vector: (Y, freq):

% Compute spectrum
Y = fft(signal);

% Perform scaling of spectrum
Y = Y/length(Y);

% Frequency vector
F_0 = fs/length(Y);
freq = 0:F_0:fs-F_0;

end