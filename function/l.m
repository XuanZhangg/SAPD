function[value] = l(i, x,A,b, mod)
%the loss function l_i(x)
%when mod == 'value', the output is the value of l_i(x)
%when mod == 'gradient', the output is the gradeint of l_i(x)
if strcmp(mod, 'value')
    if isinf(exp(-b(i)*A(i,:)*x))
        value = -b(i)*A(i,:)*x;
    else
        value = log(1 + exp(-b(i)*A(i,:)*x));
    end
elseif strcmp(mod, 'gradient')
    if isinf(exp(-b(i)*A(i,:)*x))
        value = -b(i)*A(i,:)';
    else
        value = -b(i)*exp(-b(i)*A(i,:)*x)/(1 + exp(-b(i)*A(i,:)*x))*A(i,:)';
    end
end
end

