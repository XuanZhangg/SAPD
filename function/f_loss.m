function [res] = f_loss(x,y,A,b,lambda,mu_y)
%compute the objective value 
b = reshape(b,length(b),1);
mid_value = -A*x.*b;
if max(isinf( exp(mid_value) ))
    res = y'* mid_value + lambda/2*x'*x - mu_y/2*y'*y;
else
    res = y'*log( 1 + exp(mid_value) ) + lambda/2*x'*x - mu_y/2*y'*y;
end
end