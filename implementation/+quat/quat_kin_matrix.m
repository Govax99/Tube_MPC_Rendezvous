function A = quat_kin_matrix(w)
%DESCRIPTION compute matrix used in quaternion kinematics
%
% INPUT:
%    w                   angular velocity (rad/s)
%
% OUTPUT:
%	 A                   matrix quaternion kinematics
%
    A = [0, w(3), -w(2), w(1); ...
        -w(3), 0,  w(1), w(2); ...
        w(2), -w(1), 0, w(3); ...
        -w(1), -w(2), -w(3), 0];
end

