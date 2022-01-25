function [result] = SAPD(para)
%SAPD Stochastic Acceleted Primal Dual Algorithm
%refer to https://arxiv.org/abs/2111.12743
%para initial
b = para.b;
A = para.A;
[n, d] = size(A);
K = para.maxIter;
zero_vector = para.zero_vector;
m_x = para.batchsize_x;
m_y = para.batchsize_y;
rho = para.rho;%convergence rate
Rho = para.Rho;%||y||^2<= 2*Rho/n^2
R = para.R; %diameter of x
lambda= para.lambda; %regulizer for x
mu_y = para.mu_y; %regulizer for y

%stepsize
theta = para.theta;
tau = para.tau;
sigma = para.sigma;

%data initial
x = para.x;
y = para.y;
x(:, 1) = para.x0;
y(:, 1) = para.y0;
x(:, 2) = para.x0;
y(:, 2) = para.y0;
x_bar = x;
y_bar = y;
K_N(1) = 0;
K_N(2) = 0;

for k = 2:K
    %update for y
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
    batch_3 = randperm(n, m_x);
    Gradient_x = 0;
    for i = batch_3
        Gradient_x = Gradient_x + y(i, k + 1)*l(i,x(:, k),A,b,'gradient');
    end
    P_center= x(:, k) - tau*(Gradient_x*n/m_x);
    x(:, k+1) = Proj_f(P_center, lambda*tau, R);
    K_N(k + 1) = K_N(k) + rho^(-(k - 2));
    x_bar(:, k+1) = (x_bar(:, k)*K_N(k) + x(:, k+1)*rho^(-(k-2)))/K_N(k+1);
    y_bar(:, k+1) = (y_bar(:, k)*K_N(k) + y(:, k+1)*rho^(-(k-2)))/K_N(k+1);
end

result.x = x;
result.y = y;
result.x_bar = x_bar;
result.y_bar = y_bar;

for k = 1:K
    error(k) = predict_logistic(x(:, k),A,b);
    lossvalue(k) = f_loss(x(:,k),y(:,k),A,b,lambda,mu_y);
    error_bar(k) = predict_logistic(x_bar(:, k),A,b);
    lossvalue_bar(k) = f_loss(x_bar(:,k),y_bar(:,k),A,b,lambda,mu_y);
end
result.error = error;
result.lossvalue = lossvalue;
result.error_bar = error_bar;
result.lossvalue_bar = lossvalue_bar;

end


