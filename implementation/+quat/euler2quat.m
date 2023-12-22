function q = euler2quat(e, theta)
%DESCRIPTION compute quaternion from euler axis-angle
%
% INPUT:
%    e                   euler axis vector
%    theta               euler angle (rad)
%
% OUTPUT:
%    q                   quaternion q = [q_vec; q_scalar]
%
    e = e/norm(e);
    q = [e*sin(theta/2); cos(theta/2)];
end

