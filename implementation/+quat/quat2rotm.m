function R = quat2rotm(q)
%DESCRIPTION compute rotation matrix from quaternion
%
% INPUT:
%    q                   quaternion q = [q_vec; q_scalar]
%
% OUTPUT:
%	 R                   rotation matrix
%
    w = q(4);
    x = q(1);
    y = q(2);
    z = q(3);
    
    R = [1 - 2*y^2 - 2*z^2, 2*x*y - 2*w*z, 2*x*z + 2*w*y;
         2*x*y + 2*w*z, 1 - 2*x^2 - 2*z^2, 2*y*z - 2*w*x;
         2*x*z - 2*w*y, 2*y*z + 2*w*x, 1 - 2*x^2 - 2*y^2];
end
