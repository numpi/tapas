function [m, time, y] = inv_ament(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_AMENT 

timer = tic;

% FIXME: We may want to precondition this AMEn as we do in the
% non-transpose case (@see inv_amen.m). 
        
Q = ktt_infgen(R, W, ttol);
S = inv_computeS(R, W, Q, absorbing_states, shift, ttol);
QSt = round(Q - S, ttol)';

xx = amen_block_solve({ QSt }, { pi0 }, ttol, 'nswp', 1000, 'tol_exit', tol);
m = -dot(r, xx);
y = xx;

time = toc(timer);

end

