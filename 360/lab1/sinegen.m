function y = sinegen(fsamp, fsig, nsamp)
%SINEGEN Generates a sine wave.
%   fsamp - sampling frequency in Hz
%   fsig  - sine wave frequency in Hz
%   nsamp - number of samples
tsamp = 1/fsamp;
t = 0:tsamp:(nsamp-1)*tsamp;
y = sin(2*pi*fsig*t);
end