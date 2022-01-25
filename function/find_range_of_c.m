function [res] = find_range_of_c(para)
%% %Summary of this function goes here
%find the range of c given rho
%%  %initialize the parameters.
res.c_max = -inf;
res.c_min = inf;
rho = para.rho;
mu_x = para.mu_x;
mu_y = para.mu_y;
L_xx = para.L_xx;
L_yy = para.L_yy;
L_yx = para.L_yx;
tau =(1-rho)/(mu_x*rho);
%% find a feasible c for fixed rho
    cvx_clear
    cvx_precision high
    cvx_begin sdp quiet
    variables theta sigma_inv ulpha;
    minimize( 0 );
    subject to
    [sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
        (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
        (theta/rho - 1)*L_yy, 0, sigma_inv - ulpha, -theta/rho*L_yy;...
        0, -theta/rho*L_yx, -theta/rho*L_yy, ulpha/rho] >= [1e-6*eye(4)];%[1e-8*eye(3), [0,0,0]';0,0,0,0];
    theta>=0;
    sigma_inv>=0;
    cvx_end
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        evalues = eig([sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
            (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
            (theta/rho - 1)*L_yy, 0, sigma_inv - ulpha, -theta/rho*L_yy;...
            0, -theta/rho*L_yx, -theta/rho*L_yy, ulpha/rho]);
        if evalues >= 0
            c0 = ulpha/sigma_inv;
        else
            return
        end
    end
%% find min c for fixed rho.
l = 0;
r = c0;
while r - l >= 1e-9
    c = (l + r)/2;
    cvx_clear
    cvx_precision high
    cvx_begin sdp quiet
    variables theta sigma_inv;
    minimize( theta );
    subject to
    [sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
        (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
        (theta/rho - 1)*L_yy, 0, (1 - c)*sigma_inv, -theta/rho*L_yy;...
        0, -theta/rho*L_yx, -theta/rho*L_yy, c*sigma_inv/rho] >= [1e-6*eye(4)];%[1e-8*eye(3), [0,0,0]';0,0,0,0];
    theta>=0;
    sigma_inv>=0;
    cvx_end
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        evalues = eig([sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
            (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
            (theta/rho - 1)*L_yy, 0, (1 - c)*sigma_inv, -theta/rho*L_yy;...
            0, -theta/rho*L_yx, -theta/rho*L_yy, c*sigma_inv/rho]);
        if evalues >= 0
            res.c_min = c;
            r = c;
        else
            l = c;
        end
    else
        l = c;
    end
end
%% find max c for fixed rho.
l = res.c_min;
r = 1;

while r - l >= 1e-9
    c = (l + r)/2;
    cvx_clear
    cvx_precision high
    cvx_begin sdp quiet
    variables theta sigma_inv;
    minimize( theta );
    subject to
    [sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
        (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
        (theta/rho - 1)*L_yy, 0, (1 - c)*sigma_inv, -theta/rho*L_yy;...
        0, -theta/rho*L_yx, -theta/rho*L_yy, c*sigma_inv/rho] >= [1e-6*eye(4)];%[1e-8*eye(3), [0,0,0]';0,0,0,0];
    theta>=0;
    sigma_inv>=0;
    cvx_end
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        evalues = eig([sigma_inv + mu_y - sigma_inv/rho, (theta/rho - 1)*L_yx, (theta/rho - 1)*L_yy, 0;...
            (theta/rho - 1)*L_yx, 1/tau - L_xx, 0, -theta/rho*L_yx;...
            (theta/rho - 1)*L_yy, 0, (1 - c)*sigma_inv, -theta/rho*L_yy;...
            0, -theta/rho*L_yx, -theta/rho*L_yy, c*sigma_inv/rho]);
        if evalues >= 0
            res.c_max = c;
            l = c;
        else
            r = c;
        end
    else
        l = c;
    end
end
end

