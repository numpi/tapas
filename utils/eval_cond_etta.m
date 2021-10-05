function [m, time] = eval_cond_etta(pi0, R, W, absorbing_states, ...
				conditional_indices, ...
				algorithm, debug, tol, ttol, shift, iterative_mult, ...
				use_sinc, interval_report, batch_size, max_full_size)
			
k = length(R);
n = arrayfun(@(i) size(R{i}, 1), 1 : k);
			
Q = ktt_infgen(R, W);

if ~exist('batch_size', 'var')
    batch_size = inf;
end

if ~exist('max_full_size', 'var')
    max_full_size = 16000;
end
			
% Compute the column vector a with the exit rates for the first 
% absorbing state
j = 1;
m1 = 0; m2 = 0; time = 0;

while size(conditional_indices, 1) > 0
    jmax = min(size(conditional_indices, 1), batch_size);
    
    a = Q * ktt_ej(n, conditional_indices(1, :));
    for j = 2 : jmax
        a = round(a + Q * ktt_ej(n, conditional_indices(j, :)), ttol);
    end
    conditional_indices = conditional_indices(jmax+1:end, :);

    % Compute the inverse function Qh \ a
    if strcmp(algorithm, 'ttexpsumst')
        [nm1, time1, v] = eval_inv(pi0, a, R, W, absorbing_states, algorithm, debug, ...
                                 tol, ttol, shift, iterative_mult, use_sinc, ...
                                 interval_report, [], max_full_size);
        [nm2, time2, ~] = eval_inv(v, a, R, W, absorbing_states, algorithm, debug, ...
                                 tol, ttol, shift, iterative_mult, use_sinc, ...
                                 interval_report, [], max_full_size);
    else
        [nm1, time1, v] = eval_inv(pi0, a, R, W, absorbing_states, algorithm, debug, ...
                                 tol, ttol, shift, iterative_mult, use_sinc, ...
                                 interval_report, [], max_full_size);
        [nm2, time2, ~] = eval_inv(pi0, -v, R, W, absorbing_states, algorithm, debug, ...
                                 tol, ttol, shift, iterative_mult, use_sinc, ...
                                 interval_report, [], max_full_size);
    end
    
    m1 = m1 + nm1;
    m2 = m2 + nm2;
    
    time = time + time1 + time2;
end

m = m2 / m1;
			
end
