function [v, theta] = quat2euler(q)
%DESCRIPTION compute euler axis-angle from quaternion
%
% INPUT:
%    q                   quaternion q = [q_vec; q_scalar]
%
% OUTPUT:
%	 v                   euler axis (normalized)
%    theta               euler angle
%
theta = 2*acos(q(4));
x = q(1)/sqrt(1 - q(4)*q(4));
y = q(2)/sqrt(1 - q(4)*q(4));
z = q(3)/sqrt(1 - q(4)*q(4));
v = [x y z]';
end

