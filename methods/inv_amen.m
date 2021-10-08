function [m, time, y] = inv_amen(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_AMEN 

timer = tic;

% FIXME: we need to evaluate the number of terms in a sensible way.
expn = 8;
expsums_method = 'sinc';

n = size(pi0);

max_full_size = 16000;

DA = inv_computeDA(R, W, shift);
x0 = ttexpsummldivide(DA, -r, expn, ttol, expsums_method);

Q = ktt_infgen(R, W, ttol);
S = inv_computeS(R, W, Q, absorbing_states, shift, ttol);
QS = round(Q - S, ttol);

reg = norm(QS) * ttol * 1e2;
QS = round(QS + reg * tt_eye(n), ttol);

xx = amen_block_solve({ QS }, { r }, min(tol, 1e-8), ...
    'nswp', 1000, 'tol_exit', tol, 'kickrank', 2, ...
    'x0', x0, 'max_full_size', max_full_size);
xx = round(xx, ttol);
res = norm(QS * xx - r) / norm(r);            

m = -dot(pi0, xx);

if res > tol
    error('AMEN did not converge within the prescribed number of sweeps');
end

time = toc(timer);
y = xx;

end

