function [m, time, y] = inv_ament(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_AMENT 

timer = tic;

n = arrayfun(@(i) size(R{i}, 1), 1 : length(R));

% FIXME: We may want to precondition this AMEn as we do in the
% non-transpose case (@see inv_amen.m). 
Q = ktt_infgen(R, W, ttol);
S = inv_computeS(R, W, Q, absorbing_states, shift, ttol);
QSt = round(Q - S, ttol)';

% We regularize the system to avoid problems with singular Q
reg = norm(QSt) * ttol * 1e2;
QSt = round(QSt + reg * tt_eye(n(end:-1:1)), ttol);

max_full_size = 16000;

x0 = 0 * pi0;

xx = amen_block_solve({ QSt }, { pi0 }, min(tol, 1e-8), ...
    'nswp', 1000, 'tol_exit', tol, 'kickrank', 2, ...
    'x0', x0, 'max_full_size', max_full_size);
    
res = norm(QSt * xx - pi0) / norm(pi0);

if res > tol
    error('AMEN did not converge within the prescribed number of sweeps');
end
    
m = -dot(r, xx);
y = xx;

time = toc(timer);

end

