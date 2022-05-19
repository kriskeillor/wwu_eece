% f = letters2pam2(str)
% encode a string of ASCII text into +/-1, +/-3
% using a 7-bit ASCII encoding and a (5,2) block code

function f = letters2pam2(str);          % call as Matlab function
if mod(length(str),4)
  error('Number of characters in input text string must be a multiple of 4.') 
end
f=text2bin(str);  % convert to 7-bit ASCII
g=[1 0 1 0 1;   % channel code generator matrix
   0 1 0 1 1];
f=reshape(mod(reshape(f,2,length(f)/2)'*g,2)',length(f)*5/2,1)';
f=(4*f(1:2:end)+2*f(2:2:end))-3; % convert from bits to 4-PAM

