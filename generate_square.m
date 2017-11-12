function [time_vector, signal ] = generate_square( T_0, f_s, T)
%generate_square
    % Function call:
    % generate_square( T_0, f_s, T)
    % T_0 is the time in which the signal is high
    % T is the duration of the signal
    % f_s is sampling frequency (MUST be an integer)
    time_vector_high=[0:1/f_s:T_0-1/f_s];
    signal_high=ones(1, length(time_vector_high));
    [time_vector, signal]=zero_pad(time_vector_high, signal_high, T);
end

