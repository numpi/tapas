function [DA, scl, cnd] = inv_computeDA(R, W, shift)
%INV_COMPUTEDA 

DeltapC = diagblocks(R, W, shift);
k = length(R);

DA = R;
for j = 1 : k
    DA{j} = (- DA{j} + DeltapC{j});
end

% Compute minimum and maximum eigenvalues of a Kronecker sum
[scl, cnd] = minmaxeig(DA);

end

