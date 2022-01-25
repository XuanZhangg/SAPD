function [ result ] = OGDA(para)
%OGDA Optimistic Gradient Descent Ascent Algorithm
%refer to https://arxiv.org/abs/2002.05683
%para initial
b = para.b;
A = para.A;
[n, d] = size(A);
maxIter = para.maxIter;
L_xx = para.L_xx;
L_xy = para.L_xy;
L_yx = para.L_yx;
L_yy = para.L_yy;
L = max([L_xx, L_yy, L_yx, L_xy]);
zero_vector = para.zero_vector;
m_x = para.batchsize_x;
m_y = para.batchsize_y;
Rho = para.Rho;%||y||^2<= 2*Rho/n^2
R = para.R; %diameter of x
lambda= para.lambda; %regulizer for x
mu_x = lambda;
mu_y = para.mu_y; %regulizer for y

%stepsize
theta = 1;
tau = 1/(8*L);
sigma = tau;

%data initial
y = NaN*ones(n, maxIter);
x = NaN*ones(d, maxIter);
y_bar = NaN*ones(n, maxIter);
x_bar = NaN*ones(d, maxIter);
x_bar(:, 1) = zeros(d,1);
y_bar(:, 1) = zeros(n,1);
x(:, 1) = para.x0;
y(:, 1) = para.y0;
x(:, 2) = para.x0;
y(:, 2) = para.y0;

for k= 2 : maxIter
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
    
    P_center = y(:, k) + sigma*(1 + theta)*(Gradient_y_2*n/m_y - mu_y*y(:, k)) - sigma*theta*(Gradient_y_1*n/m_y - mu_y*y(:, k - 1));
    y(:, k+1) = Proj_g(P_center, 0, Rho); %rho is diameter here, not the convergence
    
    %update for x
    m = para.batchsize_x;
    if  k == 2
        Gradient_x_1 = 0;
        batch_3 = randperm(n, m_x);
        for i = batch_3
            Gradient_x_1 = Gradient_x_1 + y(i, k+1)*l(i,x(:, k),A,b,'gradient');
        end
    else
        Gradient_x_1 = Gradient_x_2;
    end
    batch_4 = randperm(n, m_x);
    Gradient_x_2 = 0;
    for i = batch_4
        Gradient_x_2 = Gradient_x_2 + y(i, k+1)*l(i,x(:, k),A,b,'gradient');
    end
    
    P_center= x(:, k) - (1+theta)*tau*(lambda*x(:, k) + Gradient_x_2*n/m_x) + theta* tau*(lambda*x(:, k-1) + Gradient_x_1*n/m_x);
    x(:, k+1) = Proj_f(P_center, 0, R);
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
