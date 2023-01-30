% Solve LARS problem with L1 penalty on the input time series and choose the optimal solution based
% on the AIC or BIC criterion given by the param 'criterion'.

function [x, ParametersOut] = lars(y, ParametersIn)

    % Solve LARS problem with L1 penalty and get regularization path
    [x, ParametersOut] = lars_l1(y, ParametersIn);
end