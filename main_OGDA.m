clear all
close all
run('generate_data.m')
addpath('./data')
addpath('./function')
%% Run OGDA for sim_num times
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
            
            % runOGDA
            result = OGDA(para);
            
            %record the iteration sequence if needed
            %x(:, :, iter) = result.x;
            %y(:, :, iter) = result.y;
            
            %We don't record the ergodic sequence for OGDA
            %record the train loss
            train_lossvalue(iter, :) = result.lossvalue;
            train_error(iter, :) = result.error;
            
            %record the distance to optimal point if needed. Here
            %we only compute it for dataset 3 and dataset 4
            for k = 1:para.maxIter
                if data_id == 3 || data_id == 4 || data_id == 5%(we only compute distance for one data set)
                    distance(iter, k) = norm(result.x(:, k) - x_opt)^2 + norm(result.y(:, k) - y_opt);
                end
                
                %record the test error
                test_error(iter, k) = predict_logistic(result.x(:, k),A_test,b_test);
            end
        end
        
        %save the data
        name= ['.\data\', dataname{data_id}, '\OGDA',' batchsize=', num2str(para.batchsize_x), ' sim_num =',num2str(sim_num),' data_id =', num2str(data_id), '.mat'];
        if data_id == 3 || data_id == 4 || data_id == 5
            save(name, 'train_lossvalue', 'train_error','para', 'test_error', 'distance');
        else
            save(name, 'train_lossvalue', 'train_error','para', 'test_error');
        end
    end
end