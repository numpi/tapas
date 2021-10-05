function [m, time, y] = inv_full_from_tt(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_FULL_FROM_TT 

[Q, Delta] = ktt_infgen(R, W, ttol);
S = inv_computeS(R, W, Delta, absorbing_states, shift, ttol);

time = tic;

QS = round(Q - S, ttol);
xx = full(QS) \ full(r);
xx = tt_tensor(xx, ttol, size(r));
m = -dot(pi0, xx);

time = toc(time);
y = xx;

end

