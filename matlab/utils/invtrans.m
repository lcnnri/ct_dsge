function[x] = invtrans(x,trspec)
if isempty(trspec) | trspec==1
    x = x^2/(1+x^2);
elseif trspec == 2
    x = exp(x);
    a=1e-5; b=10; c=100;
%     x = a + exp(c*(x-b));
%     x=abs(x);
%     x=x;
end
