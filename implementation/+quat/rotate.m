function vp = rotate(v, q)
%ROTATE rotate vector v with quaternion q (rotation about equivalent
%euler axis of euler angle magnitude)
%
% INPUT:
%    v                   vector
%    q                   quaternion q = [q_vec; q_scalar]
%
% OUTPUT:
%	 vp                  rotated vector
%
    if size(v,1)*size(v,2) ~= 3 || size(q,1)*size(q,2) ~= 4
        error("Dimensions do not match")
    end
    if isrow(q)
        q = q';
    end
    if isrow(v)
        v = v';
    end
    t = 2*cross(q(1:3), v);
    vp = v + q(4)*t + cross(q(1:3), t);
end
