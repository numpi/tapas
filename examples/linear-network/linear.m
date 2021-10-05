%

n = 20;
k = 8;

sz = k * ones(1, n);

R = cell(1, n);
W = cell(n-1, n);

for j = 1 : n
	% R{j} = tt_matrix(gallery('jordbloc', k, 0));
    R{j} = tt_matrix(triu(rand(k),1));
end

for i = 1 : size(W, 1)
    for j = 1 : i-1
        W{i,j} = tt_matrix(eye(k));
    end
    
    W{i, i}  = tt_matrix([ rand(k-1, k); zeros(1, k) ]);
    W{i,i+1} = tt_matrix(rand(k,k));
    
    for j = i+2 : n
        W{i,j} = tt_matrix(eye(k));
    end
end

absorbing_states = k * ones(1, n);
pi0 = ktt_ej(sz, ones(1, n));
r   = ktt_ones(sz) - ktt_ej(sz, k * ones(1, n));

if n <= 7
	m = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
				 'algorithm', 'spantree', ...
				 'absorbing_states', absorbing_states);
end
			 
m = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
				 'tol', 1e-4, 'algorithm', 'amen', 'ttol', 1e-10,  ...
				 'absorbing_states', absorbing_states);