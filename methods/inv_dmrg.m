function [m, time, y] = inv_dmrg(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_DMRG 

timer = tic;
		
Q = ktt_infgen(R, W, ttol);
S = inv_computeS(R, W, Q, absorbing_states, shift, ttol);
QS = round(Q - S, ttol);

% Tolerance has been experimentally adjusted to deliver reasonable
% results. 
maxswp = 100;
xx = dmrg_solve2(QS, r, tol * 1e-2, 'nswp', maxswp);

res = norm(QS * xx - r) / norm(r);

if res > tol
    error('DMRG did not converge within the prescribed number of sweeps');
end

m = -dot(pi0, xx);
y = xx;

time = toc(timer);
end

