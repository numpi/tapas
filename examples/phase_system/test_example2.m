% function test_example(n)
%TEST_EXAMPLE 

debug = true;
ttol = 1e-8;
tol  = 1e-3;

rng('default');

% Number of components
n = 5;

% Topology of failure propagation between the components
T = createTopology(n, .5, 'starnoloops');

% We assume that the local transitions are phase independent
iso = rand(1, n);
cas = rand(1, n);
mu  = rand(1, n);
ok  = rand(1, n);
red = 3;
	
[R, W, M] = CCFModel2(n, T, iso, cas, mu, ok, red);
	
% Convert everything to TT format
for i = 1 : length(R)
	R{i} = tt_matrix(R{i});
	for j = 1 : size(W, 1)
		W{j,i} = tt_matrix(W{j,i});
	end
end

absorbing_states = [ (red+2) * ones(1,n) ; (red+1) * ones(1,n) ];

sz = [ ones(1,n) * (red+2) ];
pi0 = ktt_ej(sz, ones(1,n));

% We do not need to put to zeros the rewards on absorbing states, these are
% automatically ignored by eval_measure('inv', ...). 
r   = ktt_ones(sz);

if (red+2)^n < 1e5
	m = eval_measure('inv', pi0, r, R, W, ...
		'absorbing_states', absorbing_states, ...
		'debug', debug, ...
		'algorithm', 'spantree');
end

m = eval_measure('inv', pi0, r, R, W, ...
	'absorbing_states', absorbing_states, ...
	'debug', debug, ...
	'algorithm', 'ttexpsums2', 'ttol', 1e-8);


% end

