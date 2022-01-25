function [error] = predict_logistic(x,A,b)
%predict function given training result x and data matrix A
% label b is to compute error
predict = A*x >= 0;
predict = 2*(predict - 0.5);
b = reshape(b,length(b),1);
 error = 1 - sum(b == predict)/length(b);
end

