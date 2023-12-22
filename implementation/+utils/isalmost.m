function res = isalmost(a, b, tol)
% DESCRIPTION: assert if a is almost equal to b, considering the provided tolerance
    if nargin < 3
        tol = sqrt(eps);
    end

    res = all(abs(a - b) < tol);
end

