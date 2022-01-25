function [res] = Proj_f(x, lambda, R)
%PROJ_G argmin_{p} lambda/2 ||p||^2 + 1/2*||p - p||^2
%              s.t. ||p||^2<= R^2

p = (1/(lambda+1))*x;
norm_of_p = norm(p);

if norm_of_p<= R
    res = p;
else
    res = p/norm_of_p*R;
end

%The following uses cvx
% n = length(x);
% cvx_begin quiet
% variables p(n)
% minimize( lambda/2*p'*p + 1/2*(x-p)'*(x-p))
% subject to
% norm(p)<= R;
% cvx_end
% res = p; 
end

