function qc = conj(q)
%DESCRIPTION compute conjugate of quaternion
%
% INPUT:
%    q                   quaternion q = [q_vec; q_scalar]
%
% OUTPUT:
%	 qc                  conjugate quaternion
%
    qc = zeros(size(q));
    qc(1:3) = -q(1:3);
    qc(4) = q(4);
end
