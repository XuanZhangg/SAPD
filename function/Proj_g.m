function [res] = Proj_g(p, lambda, rho)
%PROJ_G argmin_{p} lambda/2 ||p||^2 + 1/2*||p - x||^2
%              s.t. ||p-1/n||^2<= 2rho/n^2
%                    sum(p) == 1, p>=0.

p = (1/(lambda+1))*p;
u = sort(p,'descend');
n = length(u);
r = 2*rho/n^2;

b = 0;
for k = 1:n
    b = b + u(k);
    if (b - 1)/k < u(k)
        K = k;
    else
        break
    end
end
tau = (b - u(k) - 1)/K;
res = max(p - tau, 0);

if norm(res)^2 <= r + 1/n
    return
else
    u(n+1) = 0;
    a = 0;
    b = 0;
    for k = 1:n
        a = a + u(k).*u(k);
        b = b + u(k);
        c = (r + 1/n - 1/k)/(a - b^2/k);
        if c>=0
            lambda = sqrt(c);
            upper_bound = 1/(b - k*u(k));
            if k == n
                lower_bound = 0;
            else
                lower_bound =  1/(b - k*u(k+1));
            end
            if lambda>=lower_bound && lambda <= upper_bound
                res = max(0, lambda*p - ((b*lambda) - 1)/k);
                return
            end
        end
    end   
end
end

