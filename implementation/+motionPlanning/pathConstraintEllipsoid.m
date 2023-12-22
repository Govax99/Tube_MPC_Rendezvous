function [c, ceq, cGrad, ceqGrad] = pathConstraintEllipsoid(t, x, parameters, u)
% PATHCONSTRAINTELLIPSOID This function calculates the path constraints using ellipsoidal keep-out zones.
%
% INPUTS:
% t: Time vector
% x: State vector
% parameters: Structure containing the parameters of the ellipsoid and attitude
%   - parameters.ellipse.r: Radius of the ellipsoid
%   - parameters.ellipse.xc: Center of the ellipsoid
%   - parameters.ellipse.P: Matrix related to the radius of the ellipsoid
%   - parameters.attitude: Attitude of the ellipsoid
% u: Control input
%
% OUTPUTS:
% c: Inequality constraint vector
% ceq: Equality constraint vector (empty in this case)
% cGrad: Gradient of the inequality constraint vector
% ceqGrad: Gradient of the equality constraint vector (empty in this case)
%

%Sphere definition
r = parameters.ellipse.r;
xc = parameters.ellipse.xc;
P = parameters.ellipse.P;

attitude = parameters.attitude;
nTime = size(x,2);
c = zeros(nTime,nTime);
for i = 1:length(attitude)
    if (size(attitude{i},3) == nTime)
        %debugSituation(T,C,attitude{i},x(1:3,:))
        for j = 1:nTime
            R = attitude{i}(:,:,j);
            xt = R.'*x(1:3,j);
            %c(j) = 1 - ((xt(1) - xc(1))^2/r(1)^2 + (xt(2) - xc(2))^2/r(2)^2 + (xt(3) - xc(3))^2/r(3)^2);
            c(j,j) = 1 - (xt - xc)'*P*(xt - xc); %-> P can be recreated from r
        end
        break;
    end
end
ceq = [];
if nargout == 4
    ceqGrad = [];
    for i = 1:length(attitude)
        if (size(attitude{i},3) == nTime)
            cGrad = zeros(nTime,10,nTime); %10 = [time + state + action];
            for j = 1:nTime
                R = attitude{i}(:,:,j);
                cGrad(j,2:4,j) = -2*R*P*(R.'*x(1:3,j) - xc);
            end
            break;
        end
    end
end

end

