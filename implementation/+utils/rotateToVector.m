function R = rotateToVector(v)
%ROTATETOVECTOR matrix to rotate [1,0,0] to the direction given by vector v
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
end

