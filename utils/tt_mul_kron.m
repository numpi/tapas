function w = tt_mul_kron(A, v)
%TT_MUL_KRON Compute w = (A{1} \otimes ... \otimes A{k}) * w in TT format.

if isa(A, 'tt_vector')
    % Alternative implementation
    ps = v.ps; core = v.core; r = v.r; n = v.n;
    for j = 1 : length(A)
        % Get the j-th core
        C = reshape(core(ps(j) : ps(j+1)-1), [r(j), n(j), r(j+1)]);
        C = reshape(permute(C, [2 1 3]), [n(j), r(j) * r(j+1)]);
        C = A{length(A)-j+1} * C;
        C = permute(reshape(C, [n(j) r(j) r(j+1)]), [2 1 3]);
        core(ps(j) : ps(j+1)-1) = C;
    end

    w = v;
    w.core = core;
else
    % This works for tt_matrices as well. 
    II = tt_matrix(A{1});
    for i = 2 : length(A)
        II = tkron(tt_matrix(A{i}), II);
    end
    w = II * v;
end


end

