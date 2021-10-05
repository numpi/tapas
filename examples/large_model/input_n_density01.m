function time = input_n_density01(n, method, topology, it_mult)
%
% Constructs the example decribed in the paper "Stochastic modeling and 
%   evaluation of large interdependent composed models through Kronecker 
%   algebra and Exponential sums", by G. Masetti, L. Robol, S. Chiaradonna,
%   F. Di Giandomenico, submitted.

if ~exist('it_mult', 'var')
	it_mult = false;
end

%  topology = full((eye(n,n)+sprand(n,n, 1 / n)) > 0)
if ~exist('topology', 'var')
    topology = createTopology(n, 0.2, 'starnoloops');
end

disp(topology);

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
  
  % Select the absorbing state -- state 3 in all the subsystems
  absorbing_states = 3 * ones(1, n);
  
  % Starting state of the system
  pi0 = ktt_ej(3*ones(1,n), ones(1,n));
  r   = round(ktt_ones(3*ones(1,n)) - ktt_ej(3*ones(1,n),3*ones(1,n)), 1e-8);

  % Compute the measure
  if n <= 10
    % m = computeMTTF(R, W, 1e-6, 1e-3, 'spantree');
	m = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
					 'algorithm', 'spantree', ...
					 'absorbing_states', absorbing_states);
				 pause
  end
  
  [m, time] = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
				     	   'algorithm', method, ...
						   'absorbing_states', absorbing_states, ...
						   'ttol', 1e-10, 'iterative_mult', it_mult);


end

