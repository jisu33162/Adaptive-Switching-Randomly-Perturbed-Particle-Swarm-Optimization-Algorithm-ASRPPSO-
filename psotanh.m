 function y = psotanh(x,a,b,c,d)
 y = 2*b ./ (1 + exp(-2*a .* (x-c))) - d;
% function y = psotanh(x)
% y = 2 ./ (1 + exp(-2 .* (x))) - 1;