function Y = ttexpsummldivide(A, B, n, tol, method, rmax)
%TTEXPSUMINV Compute A \ B using exponential sums. 

if ~exist('method', 'var')
    method = 'sinc';
end

if ~exist('rmax', 'var')
    rmax = inf;
end

[a,b] = expsums(n, method);
k = length(A);

% Replace A with the cell array of dense matrices
for i = 1 : k
    A{i} = full(A{i});
end

II = cell(1, k);
for i = 1 : k
    II{i} = expm(-b(1) * A{i});
end

Y = a(1) * tt_mul_kron(II, B);

for j = 2 : n
    for i = 1 : k
        II{i} = expm(-b(j) * A{i});
    end

    Y = Y + a(j) * tt_mul_kron(II, B);
    Y = round(Y, tol, rmax);
end


end

