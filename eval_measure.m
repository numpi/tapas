function [m, time] = eval_measure(fun, pi0, rewards, R, W, varargin)
%EVAL_MEASURE 

p = inputParser;

addOptional(p, 'algorithm', 'auto');
addOptional(p, 'absorbing_states', []);
addOptional(p, 'debug', true);
addOptional(p, 'tol', 1e-4);
addOptional(p, 'ttol', 1e-8);
addOptional(p, 'shift', 0);
addOptional(p, 'iterative_mult', false);
addOptional(p, 'use_sinc', true);
addOptional(p, 'interval_report', 50);
addOptional(p, 'conditional_indices', []);
addOptional(p, 'x0', []);
addOptional(p, 'moment', 0);
addOptional(p, 'batch_size', inf);
addOptional(p, 'max_full_size', 16000);

parse(p, varargin{:});

algorithm = p.Results.algorithm;
absorbing_states = p.Results.absorbing_states;
debug = p.Results.debug;
tol = p.Results.tol;
ttol = p.Results.ttol;
shift = p.Results.shift;
iterative_mult = p.Results.iterative_mult;
use_sinc = p.Results.use_sinc;
interval_report = p.Results.interval_report;
conditional_indices = p.Results.conditional_indices;
x0 = p.Results.x0;
moment = p.Results.moment;
batch_size = p.Results.batch_size;
max_full_size = p.Results.max_full_size;

switch fun
	case 'inv'
		[m, time] = eval_inv(pi0, rewards, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, max_full_size);
	case 'inv2'
		[~, time1, y] = eval_inv(pi0, rewards, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, max_full_size);
		[m, time2] = eval_inv(pi0, -y, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, max_full_size);
		time = time1 + time2;
        
    case 'momentk'
        y = rewards;
        time = 0.0;
        for j = 1 : moment
            [m, tm, y] = eval_inv(pi0, y, R, W, absorbing_states, ...
                algorithm, debug, tol, ttol, shift, ...
				iterative_mult, use_sinc, interval_report, x0, max_full_size);
            
            if j < moment
                % FIXME: We should probably check if the method is using
                % TT-vectors or not -- spantree may not be the only
                % exception. 
                if ~strcmp(algorithm, 'spantree')
                    y = round(rewards .* y, ttol);
                else
                    y = full(rewards) .* full(y);
                end
            end
            
            time = time + tm;
        end
        
        m = factorial(moment) * (-1)^(moment+1) * m;
    
    case 'tta_variance'
        [m1, time1, y] = eval_inv(pi0, rewards, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, max_full_size);
		[m2, time2] = eval_inv(pi0, -y, R, W, absorbing_states, ...
							 algorithm, debug, tol, ttol, shift, ...
							 iterative_mult, use_sinc, interval_report, x0, max_full_size);
        m = 2*m2 - m1^2;
		time = time1 + time2;
	
	case 'cond_etta'
		% Conditional expected time to absorption: we select as
		% absorbing_state to condition the first in the matrix
		% absorbing_states
		[m, time] = eval_cond_etta(pi0, R, W, absorbing_states, ...
								   conditional_indices, ...
								   algorithm, debug, tol, ttol, shift, ...
								   iterative_mult, use_sinc, interval_report, batch_size, max_full_size);
							   
		if debug
			fprintf('EVAL_MEASURE :: cond_etta :: measure = %e\n', m);
			fprintf('EVAL_MEASURE :: cond_etta :: time    = %f\n', time);
        end
		
	otherwise
		error('Unsupported measure');
end

end

