%% The data will be saved to data file
run('main_SAPD.m')
run('main_OGDA.m')
run('main_SMP.m')
%% the figure will be saved to figure file
addpath(genpath('./plot_tools'))
run('main_plot_train_error.m')
run('main_plot_distance.m')
run('main_plot_test_error.m')