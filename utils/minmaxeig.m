function [l,u] = minmaxeig(A)
%MINMAXEIG Compute minimum and maximum eigenvalues of a Kronecker sum.
%
% [L, U] = MINMAXEIG(A), where A is a cell array of matrices, computes the
%     minimum and maximum eigenvalue of 
%
%  kron(A{1}, I, ..., I) + kron(I, A{2}, I, ..., I) + kron(I, ..., I, A{k})
%
%     where I are identities of approxiate dimensions. The computation is
%     performed by computing the eigenvalues of the A{j}, which are
%     expected to be small matrices. 

l = 0.0;
u = 0.0;

for j = 1 : length(A)
    s = eig(full(A{j}));
    l = l + min(abs(s));
    u = u + max(abs(s));
end

end

