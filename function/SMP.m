function [ result ] = SMP(para)
%SMP stochastic mirror prox algorithm
%refer to https://arxiv.org/abs/0809.0815
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
L = max([L_xx+mu_x, L_yy+mu_y, L_yx, L_xy]);
zero_vector = para.zero_vector;
m_x = para.batchsize_x;
m_y = para.batchsize_y;
maxIter = para.maxIter;
Rho = para.Rho;%||y - 1/n||^2<= 2*Rho/n^2
R = para.R; %||x||^2 <= R^2

%stepsize
tau = 1/(sqrt(3)*L);

%data initial
x = para.x;
y = para.y;
x_bar = x;
y_bar = y;
%gap = zeros(1, maxIter);
x(:, 1) = x0;
y(:, 1) = y0;
x_bar(:, 1) = x0;
y_bar(:, 1) = y0;
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
    
    x_help = x(:, k) - tau*(lambda*x(:, k) + Gradient_x_1*n/m_x);
    y_help = y(:, k) + tau*(Gradient_y_1*n/m_y  - mu_y*y(:, k));
    y_help = Proj_g(y_help, 0, Rho); %Rho is diameter here, not the convergence
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
    
    x(:, k+1) = x_help - tau*(lambda*x_help + Gradient_x_2*n/m_x );
    y(:, k+1) = y_help + tau*(Gradient_y_2*n/m_y  - mu_y*y_help );
    y(:, k+1) = Proj_g(y(:, k + 1), 0, Rho); %rho is diameter here, not the convergence
    x(:, k+1) = Proj_f(x(:, k + 1), 0, R);
    x_bar(:, k+1) = (x_bar(:, k)*(k-1) + x(:, k+1))/k;
    y_bar(:, k+1) = (y_bar(:, k)*(k-1) + y(:, k+1))/k;
end

result.x = x;
result.y = y;
result.x_bar = x_bar;
result.y_bar = y_bar;

for k = 1:maxIter
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
