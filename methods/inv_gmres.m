function [m, time, y] = inv_gmres(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_GMRES 

n = size(pi0);

[QQ, Delta] = ktt_infgen(R, W, ttol);
[DA, scl, cnd] = inv_computeDA(R, W, shift);

S = inv_computeS(R, W, Delta, absorbing_states, shift, ttol);
QQ = round(QQ - S, ttol);

reg = norm(QQ) * ttol * 1e2;
QQ = round(QQ - reg * tt_eye(n), ttol);

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
    expinv(b), 1e2 * ttol, ...
    'restart', 1800, 'max_iters', 500, 'tol_exit', tol);

l{1} = l{1};

m = -dot(pi0, l{1});
time = toc(timer);
y = l{1};
        

end

