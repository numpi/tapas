%

% n = 6;

sz = 2 * ones(1, n);

R = cell(1, n);
W = cell(1, n);

for j = 1 : n
	R{j} = tt_matrix([ 0 1 ; 0 0 ]);
end

for j = 1 : n
	W{1, j} = tt_matrix([ 0 1 ; 1 0 ]);
end

absorbing_states = 2 * ones(1, n);
pi0 = ktt_ej(sz, ones(1, n));
r   = ktt_ones(sz) - ktt_ej(sz, 2 * ones(1, n));

if n <= 8
	m = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
				 'algorithm', 'spantree', ...
				 'absorbing_states', absorbing_states);
end
			 
m = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
				 'algorithm', 'ttexpsums2', ...
				 'interval_report', 1, ...
				 'absorbing_states', absorbing_states);
			 