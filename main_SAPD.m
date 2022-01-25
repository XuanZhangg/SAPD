clear all
%generate the dataset and the parameters, e.g, Lipschitz constant L_xx, and strongly convex modules mu_x
run('generate_data.m')
%% compute the min convergence rate given data set and the parameters if needed
%min_rho = Get_rho_SAPD(para);
%% set the convergence rate rho for every data set
%recall that dataname = {'Adult', 'MNIST','Dry Bean', 'Arcene', 'DrivFace'};
rho_set = {[0.9992], [0.9985], [0.9986 0.9997],[0.989 0.992], [0.9998]}; % Those rho are computed by Get_rho_SAPD(para) and then plus a tolerence
for rho = rho_set{data_id}
    % set the convergence rate rho
    para.rho = rho;
    
    %% Find the stepsize of SAPD given convergence rate rho
    %minimize the robustness bound for given rho, and yield the stepsize tau sigma theta
    range_of_c = find_range_of_c(para);
    c_min = range_of_c.c_min;
    c_max = range_of_c.c_max;
    grid_c = linspace(c_min, c_max,50);
    R = inf;
    for c = grid_c
        para.c = c;
        info =  min_Robustness_SAPD(para);
        if info.value < R
            R = info.value;
            info_of_optimize_R = info;
        end
    end
    para.tau = info_of_optimize_R.tau;
    para.sigma = info_of_optimize_R.sigma;
    para.theta = info_of_optimize_R.theta;
    
    %% use the above step size to run SAPD
    %batchsize is an array that can include many choices, the alg will
    %run them one by one
    for batchsize_x = batchsize
        for batchsize_y = [batchsize_x]
            %% Run SAPD for sim_num times
            for iter = 1:sim_num
                %report the cur iteration
                iter
                
                %set the random seed
                rng(nthprime(iter))
                
                % injected noise if needed
                % para.delta_x = 0;
                % para.delta_y = 0;
                % para.omega_x = para.delta_x*randn(d, para.maxIter);
                % para.omega_y = para.delta_y*randn(n, para.maxIter);
                para.batchsize_x = batchsize_x;
                para.batchsize_y = batchsize_y;
                
                % runSAPD
                result = SAPD(para);
                
                %record the iteration sequence if needed
                %x(:, :, iter) = result.x;
                %y(:, :, iter) = result.y;
                
                %record the train loss
                train_lossvalue(iter, :) = result.lossvalue;
                train_error(iter, :) = result.error;
                train_lossvalue_bar(iter, :) = result.lossvalue_bar;
                train_error_bar(iter, :) = result.error_bar;
                
                for k = 1:para.maxIter
                    %record the distance to optimal point if needed. Here
                    %we only compute it for dataset 3 and dataset 4
                    if data_id == 3 || data_id == 4 
                        distance(iter, k) = norm(result.x(:, k) - x_opt)^2 + norm(result.y(:, k) - y_opt);
                        distance_bar(iter, k) = norm(result.x_bar(:, k) - x_opt)^2 + norm(result.y_bar(:, k) - y_opt);
                    end
                    
                    %record the test error
                    test_error(iter, k) = predict_logistic(result.x(:, k),A_test,b_test);
                    test_error_bar(iter, k) = predict_logistic(result.x_bar(:, k),A_test,b_test);
                end
            end
            
            %save the data
            name= ['.\data\', dataname{data_id}, '\SAPD rho=',num2str(para.rho),' batchsize=', num2str(para.batchsize_x), ' sim_num =',num2str(sim_num),' data_id =', num2str(data_id), '.mat'];
            if data_id == 3 || data_id == 4 
                save(name, 'train_lossvalue', 'train_error', 'train_lossvalue_bar', 'train_error_bar','para', 'test_error', 'distance', 'test_error_bar', 'distance_bar');
            else
                save(name, 'train_lossvalue', 'train_error', 'train_lossvalue_bar', 'train_error_bar','para', 'test_error', 'test_error_bar');
            end
        end
    end
end