function y = plot_areaerrorbar(data, options)
options.alpha      = 0.2;
options.line_width = 1.5;
%options.error      = 'var';
options.x_axis = options.x_axis(:);
color = options.color;
options.color_area = color;
options.color_line = color;
%options.handle = figure(1);

% Computing the mean and standard deviation of the data matrix
data_mean = mean(data,1);
data_std  = std(data,0,1);
data_min = min(data,[],1);
data_max = max(data,[],1);
%if infeas
index_min = (data_min>0);
%index_min(1:40) = zeros(1,40);
%else
%    index_min = (data_min>0);
%end
% Type of error plot
%switch(options.error)
%    case 'std',
%error = data_std;
%    case 'sem', error = (data_std./sqrt(size(data,1)));
%    case 'var', error = (data_std.^2);
%    case 'c95',
%error = (data_std./sqrt(size(data,1))).*1.96;
%end

% Plotting the result
%figure(options.handle);
x_vector = [options.x_axis(index_min)', fliplr(options.x_axis(index_min)')];
%patch = fill(x_vector, [data_mean+error,fliplr(data_mean-error)], options.color_area);
%patch = fill(x_vector, [data_mean.*(2 + 2*error)./(1+error),fliplr(data_mean./(1+error))], options.color_area);
patch = fill(x_vector, [data_min(index_min),fliplr(data_max(index_min))], options.color_area);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', options.alpha);
set(gca, 'YScale','log');
hold 'on'
if strcmp(options.alg, 'OGDA')
    y = plot(options.x_axis, data_mean,'-o','MarkerIndices',1:length(data_mean)/10:length(data_mean), 'color', options.color_line,'LineWidth', options.line_width);
elseif strcmp(options.alg, 'SMD')
    y = plot(options.x_axis, data_mean,'-^','MarkerIndices',length(data_mean)/20:length(data_mean)/10:length(data_mean), 'color', options.color_line,'LineWidth', options.line_width);
else
    y = plot(options.x_axis, data_mean, 'color', options.color_line,'LineWidth', options.line_width);
end
%hold on;
end