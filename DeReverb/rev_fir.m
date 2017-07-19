function h = rev_fir(T60, sr, density)
% function h = rev_fir(T60, sr, density)
%
% Generates a pseudo-random FIR filter for use in synthesizing reverb.
% Note that the T60 is related to physical room parameters by:
% T60 = 0.161*V/(S*a) where V is volume in meters, S is the surface area, and
% a is the average absorbtion/Sabine of the wall material
%
% INPUT:
%    T60     - T60 time in seconds
%    sr      - sample rate
%    density - density (per sec) of non-zero filter taps [sr/10]
% 
% OUTPUT:
%    h       - calculated IR
%

if nargin < 3
    density = sr/10;
end

f_len = T60*sr;             % filter length
p = density/sr;             % prob of non-zero filter taps
h = floor((1+p)*rand(1,f_len)-p/2);

tau = T60/log(1000);        % time constant (-60dB ~ 1/1000)
tvec = (0:length(h)-1)/sr;  % time vec (sec)
env = exp(-tvec/tau);
h = h .* env;

