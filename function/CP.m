function [result] = SAPD(para)
%UNTITLED 此处显示有关此函数的摘要
%smoothed ringed loss with SPDC
%   此处显示详细说明
b = para.b;
A = para.A;
[n, d] = size(A);

K = para.maxIter;

L_xx = para.L_xx;
L_xy = para.L_xy;
L_yx = para.L_yx;
L_yy = para.L_yy;
lambda = para.lambda;
mu_x = lambda;
mu_y = para.mu_y;
L = max([L_xx, L_yy, L_yx, L_xy]);

y = NaN*ones(n, K);
x = NaN*ones(d, K);
y_bar = NaN*ones(n, K);
x_bar = NaN*ones(d, K);
x_bar(:, 1) = zeros(d,1);
y_bar(:, 1) = zeros(n,1);
x_bar(:, 2) = x_bar(:, 1);
y_bar(:, 2) = y_bar(:, 1);
x(:, 1) = para.x0;
y(:, 1) = para.y0;
x(:, 2) = para.x0;
y(:, 2) = para.y0;
K_N(1) = 0;
K_N(2) = 0;
zero_vector = para.zero_vector;
m_x = para.batchsize_x;
m_y = para.batchsize_y;

Rho = para.Rho;%||y||^2<= 2*Rho/n^2
R = para.R; %diameter of x
lambda= para.lambda; %regulizer for x
mu_y = para.mu_y; %regulizer for y

% mu = mu_y*(mu_x + L_xx)/2/L^2*(...
%     sqrt(1 + 4*mu_x*L^2/mu_y/(mu_x + L_xx)^2) - 1);
% theta = 1 - mu;
theta = 1;
if L_yy == 0
    beta = 1;
    theta = 1 + beta*(L_xx + mu_x)*mu_y/2/L_yx^2 - sqrt(beta^2*mu_y^2*(L_xx+mu_x)^2/4/L_yx^4 + beta*mu_x*mu_y/L_yx^2);
else
    for beta = linspace(1e-6,1-1e-6,1000)
        theta1 = 1 + beta*(L_xx + mu_x)*mu_y/2/L_yx^2 - sqrt(beta^2*mu_y^2*(L_xx+mu_x)^2/4/L_yx^4 + beta*mu_x*mu_y/L_yx^2);
        theta2 = 1 + (1 - beta)^2*mu_y^2/8/L_yy^2 - sqrt((1 - beta)^4*mu_y^4/64/L_yy^4 + (1 - beta)^2*mu_y^2/4/L_yy^2);
        theta = min(theta, max(theta1,theta2));
    end
end

tau = (1-theta)/theta/mu_x;
sigma = (1-theta)/theta/mu_y;
rho = theta;


for k = 2:K
    %update for y
    m = para.batchsize_y;
    if  k == 2
        Gradient_y_1 = zero_vector;
        batch_1 = randperm(n, m_y);
        for i = batch_1
            Gradient_y_1(i) = l(i, x(:, k-1),A,b, 'value') ;
        end
    else
        Gradient_y_1 = Gradient_y_2;
    end
    batch_2 = randperm(n, m_y);
    Gradient_y_2 = zero_vector;
    for i = batch_2
        Gradient_y_2(i) = l(i, x(:, k),A,b, 'value');
    end
    
    P_center = y(:, k) + sigma*(1 + theta)*(Gradient_y_2*n/m_y) - sigma*theta*(Gradient_y_1*n/m_y);
    y(:, k+1) = Proj_g(P_center, mu_y*sigma, Rho); %rho is diameter here, not the convergence
    
    %update for x
    m = para.batchsize_x;
    batch_3 = randperm(n, m_x);
    Gradient_x = 0;
    for i = batch_3
        Gradient_x = Gradient_x + y(i, k + 1)*l(i,x(:, k),A,b,'gradient');
    end
    P_center= x(:, k) - tau*(Gradient_x*n/m_x);
    x(:, k+1) = Proj_f(P_center, lambda*tau, R);
    K_N(k + 1) = K_N(k) + rho^(-(k - 2));
    x_bar(:, k+1) = (x_bar(:, k)*K_N(k) + x(:, k+1)*rho^(-(k-2)))/K_N(k+1);
    
end
result.x = x;
result.y = y;
result.x_bar = x_bar;
result.y_bar = y_bar;
for k = 1:K
    error(k) = predict_logistic(x(:, k),A,b);
    lossvalue(k) = f_loss(x(:,k),y(:,k),A,b,lambda,mu_y);
end
result.error = error;
result.lossvalue = lossvalue;

end


