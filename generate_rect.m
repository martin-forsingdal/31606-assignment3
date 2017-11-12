function [ s_rect ] = generate_rect( N0, N_high )
% Generates a square function of N0 elements centered at k=N0/2.
% The rect function has N_high elements equal to one, and all
% other elements are equal to zero. The signal is a column vector.
% Function call:
% function [ s_rect ] = generate_rect( N0, N_high )

% Initialize signal
s_rect=ones(1, N0);
% Number of zeros
N_zeros=N0-N_high;
% If N_zeros is odd, place the most zeros first
s_rect(1:ceil(N_zeros/2))=0;
s_rect(ceil(N_zeros/2)+1+N_high:end)=0;

end

