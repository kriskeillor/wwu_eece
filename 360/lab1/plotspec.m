function plotspec(x)
%PLOTSPEC Graphs the magnitude and phase of x
%   x - sample train
stem(abs(x));
plot(angle(x));
end