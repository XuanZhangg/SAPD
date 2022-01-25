function [ result ] = Get_rho_SAPD(para)
%% %Summary of this function goes here
%Input: Lipschitz constants and convergence rate.
%Ouptput: Best convergence rate.
%%
TIMES = 0;
left_point = 0;
right_point = 1-1e-6;
mu_x = para.mu_x;
mu_y = para.mu_y;
L_xx = para.L_xx;
L_yy = para.L_yy;
L_yx = para.L_yx;
for i = 1:1000
    if i == 1
        rho = right_point;
    else
        rho = (right_point + left_point) / 2;
    end
    TIMES = TIMES + 1;
    tau =(1-rho)/(para.mu_x*rho);
    cvx_clear
    cvx_precision high
    cvx_begin sdp quiet
    variables ulpha theta sigma_inv
    minimize( theta )
    subject to
    [sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
        (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
        (theta/rho - 1)*L_yy, 0, sigma_inv - ulpha, -theta/rho*L_yy;...
        0, -theta/rho*L_yx, -theta/rho*L_yy, ulpha/rho] >=  [1e-6*eye(4)];
    theta>=0;
    ulpha>=0;
    sigma_inv>=0;
    cvx_end
    if strcmp(cvx_status, 'Solved')%||strcmp(cvx_status, 'Inaccurate/Solved')
        evalues = eig( [sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
            (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
            (theta/rho - 1)*L_yy, 0, sigma_inv - ulpha, -theta/rho*L_yy;...
            0, -theta/rho*L_yx, -theta/rho*L_yy, ulpha/rho]);
        if evalues >= 0
            right_point = rho;
            result.theta = theta;
            result.rho = rho;
            result.ulpha = ulpha;
            result.tau = tau;
            result.sigma = 1/sigma_inv;
        else
            left_point = rho;
        end
        
    else
        left_point = rho;
    end
    if abs(right_point - left_point) < 1e-15
        break
    end
end
end

