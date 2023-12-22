function X = ThomasAlgorithm(A, B, C, D)
    %THOMASALGORITHM Solves a tridiagonal linear system using the Thomas Algorithm.
    %
    % INPUTS:
    % A: Sub-diagonal vector of the tridiagonal matrix.
    % B: Diagonal vector of the tridiagonal matrix.
    % C: Super-diagonal vector of the tridiagonal matrix.
    % D: Right-hand side vector of the linear system. (cell array)
    %
    % OUTPUT:
    % X: Solution vector of the linear system.
    %

    N = length(A);
    
    % Forward elimination
    for i = 2:N
        m = A(i) / B(i - 1);
        D{i} = D{i} - m * D{i-1};
        B(i) = B(i) - m * C(i - 1);
    end
    
    % Back substitution
    X = cell(1,N);
    X{N} = D{N} / B(N);
    
    for i = N-1:-1:1
        X{i} = (D{i} - C(i) * X{i+1}) / B(i);
    end
end

