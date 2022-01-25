%% Run SMP for sim_num times
clear all
close all
run('generate_data.m')
addpath('./data')
addpath('./function')
%% Run SMP for sim_num times
%batchsize is an array that can include many choices, the alg will
%run them one by one
for batchsize_x = batchsize
    for batchsize_y = [batchsize_x]
        for iter = 1:sim_num
            %report the cur iteration
            iter
            
            %set the random seed
            rng(nthprime(iter))
            para.batchsize_y = batchsize_y;
            para.batchsize_x = batchsize_x;
            
            % runSMP
            result = SMP(para);
            %x(:, :, iter) = result.x;
            %y(:, :, iter) = result.y;
            
            %record the train loss
            train_lossvalue(iter, :) = result.lossvalue;
            train_error(iter, :) = result.error;
            train_lossvalue_bar(iter, :) = result.lossvalue_bar;
            train_error_bar(iter, :) = result.error_bar;
            
            %record the distance to optimal point if needed. Here
            %we only compute it for dataset 3 and dataset 4
            for k = 1:para.maxIter
                if data_id == 3 || data_id == 4 %(we only compute distance for one data set)
                    distance(iter, k) = norm(result.x(:, k) - x_opt)^2 + norm(result.y(:, k) - y_opt);
                    distance_bar(iter, k) = norm(result.x_bar(:, k) - x_opt)^2 + norm(result.y_bar(:, k) - y_opt);
                end
                
                %record the test error
                test_error(iter, k) = predict_logistic(result.x(:, k),A_test,b_test);
                test_error_bar(iter, k) = predict_logistic(result.x_bar(:, k),A_test,b_test);
            end
        end
        
        %save the data
        name= ['.\data\', dataname{data_id}, '\SMP',' batchsize=', num2str(para.batchsize_x), ' sim_num =',num2str(sim_num),' data_id =', num2str(data_id), '.mat'];
        if data_id == 3 || data_id == 4
            save(name, 'train_lossvalue', 'train_error', 'train_lossvalue_bar', 'train_error_bar','para', 'test_error', 'distance', 'test_error_bar', 'distance_bar');
        else
            save(name, 'train_lossvalue', 'train_error', 'train_lossvalue_bar', 'train_error_bar','para', 'test_error', 'test_error_bar');
        end
    end
end
