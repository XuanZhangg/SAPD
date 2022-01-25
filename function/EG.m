function [ result ] = EG(para)
%ADP 此处显示有关此函数的摘要
%   此处显示详细说明
%para initial
A = para.A;
b = para.b;
[n, d] = size(A);
x0 = para.x0;
y0 = para.y0;
L_xx = para.L_xx;
L_xy = para.L_xy;
L_yx = para.L_yx;
L_yy = para.L_yy;
lambda = para.lambda;
mu_x = lambda;
mu_y = para.mu_y;
L = max([L_xx+mu_x, L_yy, L_yx, L_xy]);
omega_x = para.omega_x;
omega_y = para.omega_y;
omega_x_help = para.omega_x_help;
omega_y_help  = para.omega_y_help ;
zero_vector = para.zero_vector;
m_x = para.batchsize_x;
m_y = para.batchsize_y;
maxIter = para.maxIter;
R = para.R;
Rho = para.Rho;
R  = sqrt(d);
%stepsize
tau = 1/(4*L);
%data initial
x = para.x;
y = para.y;
%gap = zeros(1, maxIter);
x(:, 1) = x0;
y(:, 1) = y0;
k = 1;
%gap(1) = fun(x(:, 1), zero_y, 0, '' ) - fun(zero_x, y(:, 1), 0, '');

for k = 1:maxIter
    batch_1 = randperm(n, m_x);
    Gradient_x_1 = 0;
    for i = batch_1
        Gradient_x_1 = Gradient_x_1 + y(i, k)*l(i,x(:, k),A,b,'gradient');
    end
    batch_2 = randperm(n, m_y);
    Gradient_y_1 = zero_vector;
    for i = batch_2
        Gradient_y_1(i) = l(i, x(:, k),A,b, 'value') ;
    end
    x_help = x(:, k) - tau*(lambda*x(:, k) + Gradient_x_1*n/m_x + omega_x(:, k));
    y_help = y(:, k) + tau*(Gradient_y_1*n/m_y + omega_y(:, k) - mu_y*y(:, k));
    y_help = Proj_g(y_help, 0, Rho); %rho is diameter here, not the convergence
    x_help  = Proj_f(x_help, 0, R);
    
    batch_3 = randperm(n, m_x);
    Gradient_x_2 = 0;
    for i = batch_3
        Gradient_x_2 = Gradient_x_2 + y_help(i)*l(i,x_help,A,b,'gradient');
    end
    batch_4 = randperm(n, m_y);
    Gradient_y_2 = zero_vector;
    for i = batch_4
        Gradient_y_2(i) = l(i, x_help,A,b, 'value') ;
    end
    x(:, k+1) = x_help - tau*(lambda*x_help + Gradient_x_2*n/m_x + omega_x_help(:, k));
    y(:, k+1) = y_help + tau*(Gradient_y_2*n/m_y + omega_y_help(:, k) - mu_y*y_help );
    y(:, k+1) = Proj_g(y(:, k + 1), 0, Rho); %rho is diameter here, not the convergence
    x(:, k+1) = Proj_f(x(:, k + 1), 0, R);
    %k = k+1;
end
result.x = x;
result.y = y;

for k = 1:maxIter
    error(k) = predict_logistic(x(:, k),A,b);
    lossvalue(k) = f_loss(x(:,k),y(:,k),A,b,lambda,mu_y);
end
result.error = error;
result.lossvalue = lossvalue;

end
