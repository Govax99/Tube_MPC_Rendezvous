function q2 = quatInCone(q,alfa,dtheta)
%QUATINCONE generate a random quaternion given a starting one and a cone around it
%
% INPUT:
%    q                   initial quaternion q = [q_vec; q_scalar]
%    alfa                semi-angle of the cone (rad)
%    dtheta              max variation of the euler axis obtained from q

[v, theta] = quat.quat2euler(q);

%
a = [0 0 1]';
b = v;
if norm(cross(a,b)) > 1e-6
    w = cross(a,b);
    s = norm(w);
    c = dot(a,b);
    wM = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    R = eye(3) + wM + wM^2*(1-c)/s^2;
else
    R = eye(3);
end
%

n = sin(alfa);
beta = 2*pi*rand;
R2d = [cos(beta) -sin(beta); sin(beta) cos(beta)];

vp = [R2d*[1; 0]; 0];
%vp = [randn(2,1); 0];
vp = n*R*vp*rand;

v = v + vp;
v = v/norm(v);
theta = theta + dtheta*(2*rand-1);
q2 = quat.euler2quat(v,theta);
end

