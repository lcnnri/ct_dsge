function [x] = trans(x,trspec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isempty(trspec)| trspec==1
x   = sqrt(x/(1-x));
elseif trspec==2
    x = log(x);
    a=1e-5; b=10; c=100;
%     x = b + (1/c)*log(x-a);
%     x = abs(x);
%     x=x;

end



end

% invtrans