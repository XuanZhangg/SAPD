clear all
close all
run('generate_data.m')
%set the color choice array
color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

%read the batchsize from generate_data
batchsize = batchsize(1);

%read the max Iteratioon num as the length of x-axis
Leng = maxIterset(data_id);
hold on

%iter is the counter for figure
iter = 0;
%% OGDA
%sim_num = 1;
iter = iter + 1;
name = ['OGDA',' batchsize=', num2str(batchsize), ' sim_num =',num2str(sim_num),' data_id =', num2str(data_id), '.mat'];
load(name)
options.x_axis = [1:length(train_lossvalue)];
options.color = color{iter};
options.alg = 'OGDA';
p1 = plot_areaerrorbar(train_lossvalue, options);
labels{iter} = 'S-OGDA';
%% SMD
% iter = iter + 1;
% name = ['SMD',' batchsize=', num2str(batchsize), ' sim_num =',num2str(sim_num),' data_id =', num2str(data_id), '.mat'];
% load(name)
% options.x_axis = [1:length(train_lossvalue)];
% options.color = color{iter};
% options.alg = 'SMD';
% p2 = plot_areaerrorbar(train_lossvalue, options);
% labels{iter} = 'SMD';
%% SMP
iter = iter + 1;
name = ['SMP',' batchsize=', num2str(batchsize), ' sim_num =',num2str(sim_num),' data_id =', num2str(data_id), '.mat'];
load(name)
options.x_axis = [1:length(train_lossvalue)];
options.color = color{iter};
options.alg = 'SMP';
p3 = plot_areaerrorbar(train_lossvalue, options);
labels{iter} = 'SMP';
%% SAPD
%counter for SAPD figures
APD_count = 0;

%the rho set is set according to main_SAPD.m
rho_set = {[0.9992], [0.9985], [0.9986 0.9997],[0.989 0.992], [0.9998]};

for rho = rho_set{data_id}
    for batchsize_x = [batchsize]
        for batchsize_y = [1]
            iter = iter + 1;
            APD_count = APD_count + 1;
            para.rho = rho;
            
            %load the data 
            name = ['SAPD rho=',num2str(para.rho),' batchsize=', num2str(batchsize), ' sim_num =',num2str(sim_num),' data_id =', num2str(data_id), '.mat'];
            %name = ['SAPD rho=',num2str(para.rho),' batchsize=', num2str(batchsize_x), ' sim_num =',num2str(sim_num),' data_id =', num2str(data_id), '.mat'];
            load(name)
            options.x_axis = [1:length(train_lossvalue(1,:))];
            options.color = color{iter};
            %APD_figure{APD_count} = plot_areaerrorbar(train_lossvalue([1:9,11:12,14,16:18,20:38,40:48,50],:), options);%pick for batchsize = 1
            APD_figure{APD_count} = plot_areaerrorbar(train_lossvalue, options);
            
            %set the labels
            if  data_id == 3
                labels{iter} = ['SAPD ',num2str(APD_count)];
            else
                labels{iter} = ['SAPD '];
            end
        end
    end
end
hold off
%% Polish figures
%select the figure you want to show
% p1: OGDA, p2:SMD, p3: SMP
all_figures = [p1 p3];
for iter = 1:APD_count
    all_figures = [all_figures, APD_figure{iter}];
end
%create legend
legend(all_figures, labels,'Interpreter','tex','fontsize',innersize)

% create ylabel
ylabel({'$\mathcal{L}(x_k, p_k)$'},'Interpreter','latex','fontsize',outersize);

% create xlabel
xlabel({'Number of iterations: k'},'Interpreter','tex','fontsize', outersize);
box(gca,'on');
set(gca, 'XScale', 'linear','YScale', 'log', 'fontsize', innersize,'LineWidth',1)
grid on

%set the limit of x-axis and y-axis
%ylim([0 1])
xlim([0 Leng])

%% Save the figure
dataname{3} = 'DryBean';
figure_name_1 = ['./figure/',dataname{data_id}, '_obj_value_batchsize_', num2str(batchsize_x),'.fig'];
figure_name_2 = ['./figure/',dataname{data_id}, '_obj_value_batchsize_', num2str(batchsize_x),'.eps'];
saveas(all_figures,figure_name_1)
