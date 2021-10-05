function [m, time, y] = eval_inv(pi0, r, R, W, absorbing_states, ...
							  algorithm, debug, tol, ttol, shift, ...
							  iterative_mult, use_sinc, ...
                              interval_report, x0, max_full_size)
%EVAL_INV 

% FIXME: We can cleanup several parameters for eval_inv, that are not used
% anymore.

switch algorithm
    
    case 'full-from-tt'
        [m, time, y] = inv_full_from_tt(pi0, r, R, W, absorbing_states, shift, ttol, tol);
    	fprintf('m = %e (Full-From-TT), time = %f sec, max rank = %d\n', m, time, max(rank(y)));

	case 'amen'
		[m, time, y] = inv_amen(pi0, r, R, W, absorbing_states, shift, ttol, tol);		
		fprintf('m = %e (AMEn), time = %f sec\n', m, time);
        
    case 'ament'
        [m, time, y] = inv_ament(pi0, r, R, W, absorbing_states, shift, ttol, tol);		
		fprintf('m = %e (AMEnT), time = %f sec\n', m, time);
        
    case 'tt-regular-splitting'
        [m, time, y] = inv_tt_regular_splitting(pi0, r, R, W, absorbing_states, shift, ttol, tol);
        fprintf('m = %e (tt-regular-splitting), time = %f sec\n', m(1), time);       
        
    case 'dense-splitting'
        [m, time, y] = inv_dense_splitting(pi0, r, R, W, absorbing_states, shift, ttol, tol);
        fprintf('m = %e (dense-splitting), time = %f sec\n', m(1), time);        
		
	case 'dmrg'
		[m, time, y] = inv_dmrg(pi0, r, R, W, absorbing_states, shift, ttol, tol);
		fprintf('m = %e (dmrg), time = %f sec\n', m, time);        
	
	 case 'ttexpsumst'
        [m, time, y] = inv_tt_expsumst(pi0, r, R, W, absorbing_states, shift, ttol, tol);
        fprintf('m = %e (exp sums tt), time = %f sec\n', m, time);   
        
    case 'spantree'
        [m, time, y] = inv_spantree(pi0, r, R, W, absorbing_states, shift, ttol, tol);
        fprintf('m = %e (span tree), time = %f sec\n', m, time);
		
    case 'gmres'
        [m, time, y] = inv_gmres(pi0, r, R, W, absorbing_states, shift, ttol, tol);
        fprintf('m = %e (gmres), time = %f sec\n', m, time);
		
	otherwise
		error('Unsupported algorithm');
		
end



end

