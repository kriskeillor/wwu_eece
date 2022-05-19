% f = pam2letters2(seq)
% reconstruct string of +/-1 +/-3 into letters
% by employing (5,2) decoder and 7-bit ASCII 
% conversion.
function f = pam2letters2(seq)

if mod(length(seq),35)
  error('Number of PAM symbols in input must be a multiple of 35.') 
end

h=[1 0 1 0 0;
   0 1 0 1 0;
   1 1 0 0 1];
x=[0   0
   0   1
   1   0
   1   1];
cw=[0   0   0   0   0
    0   1   0   1   1
    1   0   1   0   1
    1   1   1   1   0];
syn=[0 0 0 0 0;
     0 0 0 0 1;
     0 0 0 1 0;
     0 1 0 0 0;
     0 0 1 0 0;
     1 0 0 0 0;
     1 1 0 0 0;
     1 0 0 1 0];
seq=seq(:)';
y=[(sign(seq)+1)/2; mod((seq+3)/2,2)];
y3=reshape(y(:),5,length(seq)*2/5)';
eh=mod(y3*h',2);
y4=mod(y3-syn(eh(:,1)*4+eh(:,2)*2+eh(:,3)+1,:),2);
for i=1:size(y3,1)
  y2=y4(i,:);
  for j=1:4            % recover message from codewords
    if y2==cw(j,:), z(2*i-1:2*i)=x(j,:); end
  end
end
f=bin2text(z);

