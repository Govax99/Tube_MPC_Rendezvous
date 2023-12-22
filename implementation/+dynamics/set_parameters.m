function [J_C,J_T,m_C,OM,pberth,kP_tr,kD_tr,kP_rot,u_lim,r2] = set_parameters(parameters)
%DESCRIPTION extract parameters from structure (if they exists otherwise set defaults)
%
% INPUT:
%    parameters          structure containing parameters for the dynamics
%
% OUTPUT:
%	 J_C  	             moment of inertia of chaser (in C-frame)
%	 J_T  	             moment of inertia of target (in T-frame)
%	 m_C  	             total mass of chaser
%    OM                  orbital circular velocity
%    pberth              position of berthing point in T-frame
%    kP_tr               proportional coefficient matrix (PD translation control)
%    kD_tr               differential coefficient matrix (PD translation control)
%    kP_rot              proportional coefficient matrix (P rotational control)
%    u_lim               control action limits (force limits in L-frame, torque limits in C-frame)
%    r2                  radius squared of the keep-out sphere
%
% LEGEND FRAMES:
%    C-frame             relative frame, attached to chaser object, center on chaser cm
%    T-frame             relative frame, attached to target object, center on target cm
%
if (~isfield(parameters,'J_C'))
    J_C = eye(3);
else
    J_C = parameters.J_C;
end

if (~isfield(parameters,'J_T'))
    J_T = eye(3);
else
    J_T = parameters.J_T;
end

if (~isfield(parameters,'m_C'))
    m_C = 1;
else
    m_C = parameters.m_C;
end

if (~isfield(parameters,'OM'))
    OM = 0.005;
else
    OM = parameters.OM;
end

if (~isfield(parameters,'pberth'))
    pberth = [5 0 0]';
else
    pberth = parameters.pberth;
end

if (~isfield(parameters,'kP_tr'))
    kP_tr = 0.1*eye(3);
else
    kP_tr = parameters.kP_tr;
end

if (~isfield(parameters,'kD_tr'))
    kD_tr = 1*eye(3);
else
    kD_tr = parameters.kD_tr;
end

if (~isfield(parameters,'kP_rot'))
    kP_rot = 1*eye(3);
else
    kP_rot = parameters.kP_rot;
end

if (~isfield(parameters,'u_lim'))
    u_lim = ones(1,6);
else
    u_lim = parameters.u_lim;
end

if (~isfield(parameters,'r2'))
    r2 = 3^2;
else
    r2 = parameters.r2;
end

end