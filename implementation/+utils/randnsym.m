function X = randnsym(N)
%RANDNSYM generate random symmetric matrix
X = randn(N);
for i = 1:N
    for j = 1:N
        if i > j
            X(i,j) = X(j,i);
        end
    end
end
end

