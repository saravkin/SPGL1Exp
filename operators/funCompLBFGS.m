function [f g] = funCompLBFGS(x, b, funForward, funPenalty, params)
    r = b - funForward(x, [], params);
    [f v] = funPenalty(r, params);
    g = funForward(x, -v, params);
end


% Old call: 
% spgl1(@NLfunForward, b, TAU,...) 

% New call: 
%  funObj = @(x)funCompLBFGS(x, b, @NLfunForward, @funLS, params);
%  X = lbfgs(funObj);
