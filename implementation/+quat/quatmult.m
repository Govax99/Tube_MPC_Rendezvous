function q3 = quatmult(q1, q2)
%QUATMULT compute multiplication of quaternions q1 and q2
%
% INPUT:
%    q1                  first quaternion q = [q_vec; q_scalar]
%    q2                  second quaternion q = [q_vec; q_scalar]
%
% OUTPUT:
%	 q3                  quaternion, product of q1, q2
%
    q3 = [ q1(4)*q2(1:3) + q2(4)*q1(1:3) + cross(q1(1:3),q2(1:3)); ...
            q1(4)*q2(4) - q1(1:3)'*q2(1:3); ...
        ];
end

