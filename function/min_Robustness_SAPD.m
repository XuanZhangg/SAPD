function [res] = min_Robustness_SAPD(para)
%%
%MIN_ROBUSTNESS_SAPD Summary of this function goes here:
%Compute the min robusentess.

%Detailed explanation goes here
%The inputs are Lipschitz constants and convergence rate;
%The outputs are a structure including the min robustness and stepsizes.
%Given rho and c, compute the min robustness
%%  %initialize the parameters.
res.value = inf;
rho = para.rho;
mu_x = para.mu_x;
mu_y = para.mu_y;
L_xx = para.L_xx;
L_yy = para.L_yy;
L_yx = para.L_yx;
L_xy = para.L_xy;
delta_x = para.delta_x;
delta_y = para.delta_y;
delta = para.delta;
tau =(1-rho)/(mu_x*rho);
c = para.c;
%% find min theta for fixed rho.
theta_min = NaN;
cvx_clear
cvx_precision high
cvx_begin sdp quiet
variables theta sigma_inv;
minimize( theta );
subject to
[sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
    (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
    (theta/rho - 1)*L_yy, 0, (1 - c)*sigma_inv, -theta/rho*L_yy;...
    0, -theta/rho*L_yx, -theta/rho*L_yy, c*sigma_inv/rho] >= [1e-8*eye(4)]%, [0,0,0]';0,0,0,0];
theta>=0;
sigma_inv>=0;
cvx_end
if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    evalues = eig([sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
        (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
        (theta/rho - 1)*L_yy, 0, (1 - c)*sigma_inv, -theta/rho*L_yy;...
        0, -theta/rho*L_yx, -theta/rho*L_yy, c*sigma_inv/rho]);
    if evalues >= 0
        theta_min = theta;
    end
end
%% find max theta for fixed rho.
theta_max = theta_min;
cvx_clear
cvx_precision high
cvx_begin sdp quiet
variables theta sigma_inv;
maximize( theta );
subject to
[sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
    (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
    (theta/rho - 1)*L_yy, 0, (1 - c)*sigma_inv, -theta/rho*L_yy;...
    0, -theta/rho*L_yx, -theta/rho*L_yy, c*sigma_inv/rho] >= [1e-8*eye(4)];
theta>=0;
sigma_inv>=0;
cvx_end
if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    evalues = eig([sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
        (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
        (theta/rho - 1)*L_yy, 0, (1 - c)*sigma_inv, -theta/rho*L_yy;...
        0, -theta/rho*L_yx, -theta/rho*L_yy, c*sigma_inv/rho]);
    if evalues >= 0
        theta_max = theta;
    end
end
%theta_max = min(theta_max, rho);
%%  For each theta, find the min sigma and compute the robustness.
for theta_iter = linspace(theta_min, theta_max, 5000)
    l = 0;
    r = 1/sigma_inv;
    err = r - l;
    while err > 1e-6
        sigma = (l + r)/2;
        evalues = eig([sigma^(-1) + mu_y - sigma^(-1)/rho, (theta_iter/rho - 1)*L_yx, (theta_iter/rho - 1)*L_yy, 0;...
            (theta_iter/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta_iter/rho*L_yx;...
            (theta_iter/rho - 1)*L_yy, 0, (1 - c)*sigma^(-1), -theta_iter/rho*L_yy;...
            0, -theta_iter/rho*L_yx, -theta_iter/rho*L_yy, c*sigma^(-1)/rho]);
        
        if evalues >= 0
            r = sigma;
            ulpha = c*1/sigma;
            
            %we have found ulpha and sigma, then compute the robustness
            cur_res = 2*rho/(1 - rho) * max(tau, sigma/(1 - ulpha*sigma))*1/delta^2*...
                (...
                tau/(1 + tau*mu_x)*(1/2 + 1 + sigma*theta_iter*(1+theta_iter)*L_yx/2/(1 + sigma*mu_y))*delta_x^2 ...
                +...
                sigma/(1 + sigma*mu_y)*...
                (...
                1/2 + theta_iter +...
                (1 + 2*theta_iter + (theta_iter + sigma*theta_iter*(1+theta_iter)*L_yy)/(1 + sigma*mu_y) + tau*sigma*theta_iter*(1 + theta_iter)*L_yx*L_xy/(1 + tau*mu_x)/(1 + sigma*mu_y))...
                *(1 + 2*theta_iter) ...
                + tau*theta_iter*(1 + theta_iter)*L_yx/2/(1 + tau*mu_x)...
                )*delta_y^2 ...
                );
            if cur_res < res.value
                res.value = cur_res;
                res.sigma = sigma;
                res.ulpha = ulpha;
                res.tau = tau;
                res.theta = theta_iter;
                res.rho = rho;
            end
        else
            l = sigma;
        end
        err = r - l;
    end
end
end




