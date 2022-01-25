clear all
addpath(genpath('./data'))
addpath(genpath('./function'))
dataname = {'Adult', 'MNIST','DryBean', 'Arcene', 'DrivFace'};
%% Choose the data set
%1:Adult data set, 2:MNist, 3: Dry bean, 4: Arcene,5:DrivFace
data_id = 4;
%% Set the fontsize for figure
%inner size is for legends
%outersize is for labels
innersizeset = [16,16,16,16,16];
outersizeset = [18,18,18,18,18];
innersize = innersizeset(data_id);
outersize = outersizeset(data_id);
%% load the dataset
% The comment part is to process the raw data set. The
% '' .mat'' files are processed data.
if data_id == 1
    % Adult data set
    load('Adult.mat')
%     data = readmatrix('numerical_adult_onehotcode.csv');
%     b = data(:,end);
%     A = data(:,1:end-1);
%     b(b==0)=-1;
%     [n, d] = size(A);
%     for i = [1 26 9 61 62 63]%nomrlize dense data
%         %A(:, i) = A(:, i)/norm(A(:, i));
%         A(:, i) = (A(:, i)-min(A(:, i) )) /(max(A(:, i) ) - min(A(:, i) ));
%     end
%     cv = cvpartition(size(A,1),'HoldOut',0.3);
%     idx = cv.test;
%     %Separate to training and test data
%     A_test = A(idx,:);
%     b_test = b(idx,:);
%     A = A(~idx,:);
%     b = b(~idx,:);
%     [n, d] = size(A);
elseif data_id == 2
    % MINIST DIGIT DATA SET
    digit = 4;
    load('MINIST_4.mat')
        [data, labels] = readMNIST('train-images.idx3-ubyte', 'train-labels.idx1-ubyte', 10000, 0);
        [h, l, n] = size(data);
        d =h*l;
        for i  = 1:n
            A(i, :) = reshape(data(:, :, i), 1, d);
            if labels(i) == digit
                b(i) = 1;
            else
                b(i) = -1;
            end
        end
    
        [data, labels] = readMNIST('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte', 10000, 0);
        [h, l, n] = size(data);
        d =h*l;
        for i  = 1:n
            A_test(i, :) = reshape(data(:, :, i), 1, d);
            if labels(i) == digit
                b_test(i) = 1;
            else
                b_test(i) = -1;
            end
        end
elseif data_id == 3
    load('Dry_bean.mat')
    %     data = readmatrix('Kinds_of_Seed.csv');
    %     b = data(:,end);
    %     A = data(:,1:end-1);
    %     [n, d] = size(A);
    %     for i = 1:d
    %         A(:, i) = (A(:, i)-min(A(:, i) )) /(max(A(:, i) ) - min(A(:, i) ));
    %     end
    %     cv = cvpartition(size(A,1),'HoldOut',0.3);
    %     idx = cv.test;
    %     % Separate to training and test data
    %     A_test = A(idx,:);
    %     b_test = b(idx);
    %     %A = A./max(A);
    %     %A = normalize(A(~idx,:));
    %     A = A(~idx, :);
    %     b = b(~idx);
    %     [n, d] = size(A);
elseif data_id == 4
    load('Arcene.mat')
    data = readmatrix('arcene_valid.csv');
    A_test = data(:,1:end-1);
    b_test = data(:,end);
    [row, col] = find(isnan(A_test));
    A_test(row, :) = [];
    b_test(row) = [];
    [n1, d] = size(A_test);
    data = readmatrix('arcene_train.csv');
    A= data(:,1:end-1);
        [row, col] = find(isnan(A));
    b = data(:,end);
     A(row, :) = [];
    b(row) = [];
    [n, d] = size(A);
        for i = 1:d
            if max(A_test(:, i)) == min(A_test(:, i))
                A_test(:, i) = zeros(n1,1);
            else
                A_test(:, i) = (A_test(:, i) - min(A_test(:, i)))/(max(A_test(:, i)) - min(A_test(:, i)));
            end
            if max(A(:, i)) == min(A(:, i))
                A(:, i) = zeros(n,1);
            else
                A(:, i) = (A(:, i) - min(A(:, i)))/(max(A(:, i)) - min(A(:, i)));
            end
        end
     A = A/sqrt(n);
     A_test = A_test/sqrt(n);
elseif data_id == 5
    load('drivFace_2.mat');
%     [n, d] = size(A);
%     for i = 1:d
%         A(:, i) = (A(:, i)-min(A(:, i) )) /(max(A(:, i) ) - min(A(:, i) ));
%     end
%     b(b==1) = 1;
%     b(b~=1)=-1;
%     %A = A/norm(A);
%     cv = cvpartition(size(A,1),'HoldOut',0.3);
%     idx = cv.test;
%     % Separate to training and test data
%     A_test = A(idx,:);
%     b_test = b(idx);
%     %A = A./max(A);
%     %A = normalize(A(~idx,:));
%     A = A(~idx, :);
%     b = b(~idx);
%     [n, d] = size(A);
end
%% set running time, batchsize, and simulation num
maxIterset = [3000 3000 20000 10000 10000];
para.maxIter = maxIterset(data_id);
%it is a set. One can set many choices and the algoirthm will run them one by one 
batchsize = [10];
sim_num = 1;
%% set regulizer and strongly convex parameter
if data_id == 1 %Seed %adult mu_y = 1, rho =0.9997+, mu_y = 10, rho = 0.9992+
    para.lambda = 0.01;%(good choice: 0.01)
    para.mu_x = 0.01;
    para.mu_y = 10;%(good choice: 10)
elseif data_id == 2 %MINIST
    para.lambda = 0.1;%(good choice: 0.1)
    para.mu_x = 0.1;
    para.mu_y = 10;%(good choice: 10)
elseif data_id == 3
    para.lambda = 0.01;%(good choice: 0.01)
    para.mu_x = 0.01;
    para.mu_y = 10;%(good choice: 10)
elseif data_id == 4
    para.lambda = 0.02;%(good choice: 1e-6)
    para.mu_x = 0.02;
    para.mu_y = 10;%(good choice: 10)
elseif data_id == 5
    para.lambda = 0.01;%(good choice: 1e-6)
    para.mu_x = 0.01;
    para.mu_y = 10;%(good choice: 10)
end
%% set projection diameter
%||y - 1/n||^2<= 2*Rho/n^2
para.Rho = sqrt(n);
%||x||^2 <= R^2
para.R = sqrt(d);
%% Set the variance of the noise
% In fact, if one do not inject noise, the value of delta won't have an
% effect on the result. However, one must give values to them and the
% should be equal to each other.
para.delta_y = 10;
para.delta_x = 10;
para.delta = 10;
%% Find lipshitz constant and other parameters
%minimize the robustness bound for given rho, and yield the stepsize tau sigma theta
%the function is lambda/2*||x||^2 + sum_{i=1 to n} y_i*l_i(x) -
%mu_y/2*||y||^2
%l_i is logistic loss in the codes
para.A = A;
para.b = b;
para.L_yx = norm(A);
para.L_xy = norm(A);
para.L_yy = 0;
para.L_xx = 1/4*max(sqrt(sum(A.^2, 2)));
%% Initialize iteration sequence
% d: dimension of x, n: dimension of y
para.x = zeros(d, para.maxIter);
para.y = zeros(n, para.maxIter);
para.x0 = 2*ones(d, 1);
para.y0 = 1/n*ones(n, 1);
%this is just a helper zero-vector
para.zero_vector = zeros(n, 1);