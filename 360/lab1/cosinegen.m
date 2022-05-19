function y = cosinegen(fsamp, fsig, nsamp)
%SINEGEN Generates a cosine wave.
%   fsamp - sampling frequency
%   fsig  - sine wave frequency
%   nsamp - number of samples
tsamp = 1/fsamp;                % sampling frequency
t = 0:tsamp:(nsamp-1)*tsamp;    % sampling times 
y = cos(2*pi*fsig*t);           % generated sine wave
end