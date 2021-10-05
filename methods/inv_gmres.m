function [m, time, y] = inv_gmres(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_GMRES 

k = length(R);

QQ = ktt_infgen(R, W, ttol);
[DA, scl, cnd] = inv_computeDA(R, W, shift);

b = r;

if scl ~= 0
    cnd = cnd / scl;
else
    cnd = 0;
    scl = 1;
end

QQ = QQ / scl;
for j = 1 : length(DA)
    DA{j} = -DA{j} / scl;
end
b = b / scl;

% expinv = @(x) tt_tensor(reshape(full(-kronSum(DA{:}, ttol)) \ full(x), n));
expinv = @(x) -ttexpsummldivide(DA, x, 4, ttol);
%DeltapB = round(Delta - kronSum(DeltapC{:}, ttol), ttol);
%expinv = @(x) gmres_preconditioner(DA, DeltapB, x);

% A good preconditioner ?
% M = round(kronSum(R{:}, ttol) - Delta, ttol) / scl;
% MM = round(M' * M, tol);
% expinv = @(x) amen_solve2(MM, M' * x, 1e-2);

timer = tic;
l = tt_gmres_block(...
    @(x,ttol) { round(expinv(QQ*x{1}), ttol) }, ...
    expinv(b), 1e-6, ...
    'restart', 1800, 'max_iters', 500, 'tol_exit', tol);
m = -dot(pi0, l{1});
time = toc(timer);
y = l{1};
        

end

