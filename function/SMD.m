function [ result ] = SMD(para)
%SMD stochastic mirror descent algorithm
%refer to https://www.semanticscholar.org/paper/Robust-Stochastic-Approximation-Approach-to-Nemirovski-Juditsky/96167ed3ebc9a2c3270f6ae96043e6f086eed4de
%para initial
b = para.b;
A = para.A;
[n, d] = size(A);
K = para.maxIter;
zero_vector = para.zero_vector;
Rho = para.Rho;%||y - 1/n||^2<= 2*Rho/n^2
R = para.R; %||x||^2 <= R^2
lambda= para.lambda; %regulizer for x
mu_y =para.mu_y;

%stepsize
max_a = sqrt(max(sum(A.^2, 2)));
%we used ||x||^2<R, ||y||^2<(1/n + 2*Rho/n^2), log(1+exp(x))<= log(2) + x,
%(a+b)^2<=2(a^2+b^2)
M = 2*R^2*...
   (lambda^2*R^2 + n^2/para.batchsize_x*(1/n + 2*Rho/n^2)*max_a)...
    + 2*(1/n + 2*Rho/n^2) *...
    2*(mu_y^2*(1/n + 2*Rho/n^2) +  n^2/para.batchsize_y*(log(2) + R* max_a)^2 );
gama = 2/sqrt(M*5*K);

%data initial
y = NaN*ones(n, K);
x = NaN*ones(d, K);
y_bar = NaN*ones(n, K);
x_bar = NaN*ones(d, K);
x_bar(:, 1) = para.x0;
y_bar(:, 1) = para.y0;
x(:, 1) = para.x0;
y(:, 1) = para.y0;

for k= 1 : K
    %update for y
    m_y = para.batchsize_y;
    batch_1 = randperm(n, m_y);
    Gradient_y = zero_vector;
    for i = batch_1
        Gradient_y(i) = l(i, x(:, k),A,b, 'value');
    end
    P_center = y(:, k) + gama*(Gradient_y*n/m_y - mu_y*y(:, k));
    %M_Y = max(M_Y,norm(Gradient_y*n/m_y - mu_y*y(:, k))^2);
    y(:, k+1) = Proj_g(P_center, 0, Rho); %rho is diameter here, not the convergence
    
    %update for x
    m_x = para.batchsize_x;
    batch_2 = randperm(n, m_x);
    Gradient_x = 0;
    for i = batch_2
        Gradient_x = Gradient_x + y(i, k)*l(i,x(:, k),A,b,'gradient');
    end
    
    P_center= x(:, k) - gama*(lambda*x(:, k) + Gradient_x*n/m_x);
    %M_X = max(M_X,norm(lambda*x(:, k) + Gradient_x*n/m_x)^2);
    x(:, k+1) = Proj_f(P_center, 0, R);
    
    x_bar(:, k+1) = (x_bar(:, k)*(k-1) + x(:, k+1))/k;
    y_bar(:, k+1) = (y_bar(:, k)*(k-1) + y(:, k+1))/k;
end

result.x = x;
result.y = y;
result.x_bar = x_bar;
result.y_bar = y_bar;
%M_X = 2*R^2*M_X;
%M_Y = 2*(1/n + 2*Rho/n^2) *M_Y;

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
%result.M_Y=M_Y;
%result.M_X=M_X;
end
