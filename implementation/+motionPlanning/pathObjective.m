function [obj, objGrad] = pathObjective(u)
% PATHOBJECTIVE This function computes the objective function and its
% gradients for the trajectory generation
%
% INPUTS:
% u: Control input vector
%
% OUTPUTS:
% obj: Objective function value
% objGrad: Gradient of the objective function
%

obj = 1/2*dot(u,u,1);

if nargout == 2  % Analytic gradients
    nTime = size(u,2);
    
    objGrad = zeros(10,nTime); %10 = [time + state + action];
    
    objGrad(8:10,:) = u;  %gradient obj wrt u
    
end

end

