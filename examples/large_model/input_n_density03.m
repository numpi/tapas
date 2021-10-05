function time = input_n_density02(n, method, topology)
%

if ~exist('topology', 'var')
    topology = createTopology(n, 0.5, 'starnoloops');
end

disp(topology)

  pp = symrcm(topology);

  lambda_min = 0.5;
  lambda_max = 1.5;
  lambda = lambda_min + (lambda_max-lambda_min)*rand(n,1)

  mu_min = 2000;
  mu_max = 3000;
  mu = mu_min + (mu_max-mu_min)*rand(n,1)

  p_min = 0.95;
  p_max = 1.0;
  p = p_min + (p_max-p_min)*rand(n,1)

  % Construct Rs and Ws
  [R, W] = largeModel(n, topology(pp,pp), lambda(pp), mu(pp), p(pp));
  
  % Add another absorbing state
  for j = 1 : n
	W{end+1,j} = zeros(3,3); W{end,j}(2,1) = -mu(pp(j));
	W{end,j} = tt_matrix(W{end,j});
	for i = 1 : n
		if i ~= j
			W{end,i} = zeros(3,3); W{end,i}(2,2) = 1;
			W{end,i} = tt_matrix(W{end,i});
		end
	end
  end
    
  % Select the absorbing state -- state 3 in all the subsystems
  absorbing_states = [ 3 * ones(1, n) ; 2 * ones(1, n) ];
  
  % Starting state of the system
  pi0 = ktt_ej(3*ones(1,n), ones(1,n));
  r   = round(ktt_ones(3*ones(1,n)) - ktt_ej(3*ones(1,n),3*ones(1,n)) - ktt_ej(3*ones(1,n),2*ones(1,n)), 1e-8);

  % Compute the measure
  if n <= 9
    % m = computeMTTF(R, W, 1e-6, 1e-3, 'spantree');
	m = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
					 'algorithm', 'spantree', ...
					 'absorbing_states', absorbing_states);
  end
  
  [m, time] = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
				     	   'algorithm', method, ...
						   'absorbing_states', absorbing_states, 'tol', 1e-8, ...
						   'ttol', 1e-8, 'shift', 1);


end

