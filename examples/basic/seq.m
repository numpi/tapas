%

% n = 6;

sz = 2 * ones(1, n);

R = cell(1, n);
W = cell(n, n);

for j = 1 : n
	R{j} = tt_matrix([ 0 0 ; 0 0 ]);
end

a = 0;

for i = 1 : n
	for j = 1 : i - 1
		W{i,j} = tt_matrix([ a 0 ; 0 1 ]);
	end
	
	if i > 1
		W{i,i-1} = tt_matrix([ 0 0 ; 0 1 ]);
	end
	
	W{i,i} = tt_matrix([ 0 1 ; 0 0 ]);
	
	for j = i+1 : n
		W{i, j} = tt_matrix([ 1 0 ; 0 1 ]);
	end
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
				 'algorithm', 'tt-regular-splitting', 'ttol', 1e-12,  ...
				 'absorbing_states', absorbing_states);
             
m = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
    'algorithm', 'gmres', 'ttol', 1e-12, 'tol', 1e-4, ...
    'absorbing_states', absorbing_states);