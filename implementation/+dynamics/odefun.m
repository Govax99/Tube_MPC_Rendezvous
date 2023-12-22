function [dx, dxGrad] = odefun(t, x, parameters, u)
%DESCRIPTION state-space form of chaser and target dynamics (can handle control input
% as matlab function, useful for PID control)
% (linear dynamics in L-frame, angular dynamics in C-frame, T-frame)
%
% INPUT:
%    time                current integration time
%    state               current total state [p_LC_L, v_LC_L, q_LC, w_IC_C, q_LT, w_IT_T]
%    parameters          structure containing parameters for the dynamics
%    control             control action [f, tau]
%
% OUTPUT:
%	 dx  	             derivative of state
%    f                   control force (in L-frame)
%    tau                 control torque (in C-frame)
%
% LEGEND FRAMES:
%    L-frame             LVLH, local vertical-local horizontal frame, center on target cm
%    C-frame             relative frame, attached to chaser object, center on chaser cm
%    T-frame             relative frame, attached to target object, center on target cm
%
% MONOGRAM NOTATION:
%    p_XY_Z              p: symbol of physical quantity,
%                        X: "measured from",
%                        Y: target point/frame,
%                        Z: "expressed in" frame
[~,~,m_C,OM] = dynamics.set_parameters(parameters);

p_LC_L = x(1:3,:);
v_LC_L = x(4:6,:);

% manage control
if (nargin < 4)
    try
        f = zeros(3,size(u,2));
    catch
        f = zeros(3,1);
    end
else
    f = u(1:3,:);
end

dx = [v_LC_L; ...
    (2*OM*v_LC_L(2,:) + 3*OM^2*p_LC_L(1,:) + f(1,:)/m_C); ...
    (-2*OM*v_LC_L(1,:) + f(2,:)/m_C); ...
    (-OM^2*p_LC_L(3,:) + f(3,:)/m_C); ...
    ];

if nargout == 2   % Analytic gradients
    nTime = size(u,2);
    
    %gradient dp wrt d(all)
    dpxGrad = zeros(1,10,nTime); %10 = [time + state + action];
    dpxGrad(1,5,:) = 1;
    dpyGrad = zeros(1,10,nTime);
    dpyGrad(1,6,:) = 1;
    dpzGrad = zeros(1,10,nTime);
    dpzGrad(1,7,:) = 1;
    
    %gradient dv wrt d(all)
    dvxGrad = zeros(1,10,nTime);  %10 = [time + angle + rate + torque] for every vi;
    dvxGrad(1,2,:) = 3*OM^2;
    dvxGrad(1,6,:) = 2*OM;
    dvxGrad(1,8,:) = 1/m_C;

    dvyGrad = zeros(1,10,nTime);
    dvyGrad(1,5,:) = -2*OM;
    dvyGrad(1,9,:) = 1/m_C;

    dvzGrad = zeros(1,10,nTime);
    dvzGrad(1,4,:) = -OM^2;
    dvzGrad(1,10,:) = 1/m_C;
    
    dxGrad = cat(1, dpxGrad, dpyGrad, dpzGrad, dvxGrad, dvyGrad, dvzGrad);
    
end

end

