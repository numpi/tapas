% function test_example(n)
%TEST_EXAMPLE 

debug = true;
ttol = 1e-8;
tol  = 1e-3;

rng('default');

% Number of components
n = 12;

% Number of phases in the system
nphases = 3;

% Topology of failure propagation between the components
T = createTopology(n, 1/nphases, 'rand');

% Topology and rates of the phases
TP = [ 0 1 0 ; ...
	   0 0 1 ; ...
	   1 0 0   ];
   
% We assume that the local transitions are phase independent
iso = rand(1, n);

WW = cell(1, nphases);

for j = 1 : nphases
	cas = rand(1, n);
	mu  = zeros(1, n);
	
	[R, WW{j}, M] = CCFModel(n, T, iso, cas, mu);
	% Q = infgen(R, W, M, 1e-8, false);
end

% We build the R and W factors of the large system, including the phases 
% as the first component. 
R = { TP, R{:} };

nsyncs = 0;
for j = 1 : length(WW)
	nsyncs = nsyncs + size(WW{j}, 1);
end

W = cell(1 + nsyncs, n + 1);
counter = 0;
for j = 1 : length(WW)
	S = zeros(nphases); S(j,j) = 1;
	for i = 1 : size(WW{j}, 1)
		W{i+counter, 1} = S;
		for k = 1 : n
			W{i+counter, k+1} = WW{j}{i,k};
		end
	end
	counter = counter + size(WW{j}, 1);
end

W{end,1} = -TP;
for j = 1 : n
	W{end,j+1} = [ 0 0 ; ...
		           0 1 ];
end

% Convert everything to TT format
for i = 1 : length(R)
	R{i} = tt_matrix(R{i});
	for j = 1 : size(W, 1)
		W{j,i} = tt_matrix(W{j,i});
	end
end

% Q = infgen(R, W, M, 1e-8, false);

absorbing_states = [ 1 : nphases ; 2 * ones(n, nphases) ]';

sz = [ nphases, ones(1,n) * 2 ];
pi0 = ktt_ej(sz, ones(1,n+1));

% We do not need to put to zeros the rewards on absorbing states, these are
% automatically ignored by eval_measure('inv', ...). 
r   = ktt_ones(sz);

% m = eval_measure('inv', pi0, r, R, W, ...
% 	'absorbing_states', absorbing_states, ...
% 	'debug', debug, ...
% 	'algorithm', 'spantree');

m = eval_measure('inv', pi0, r, R, W, ...
	'absorbing_states', absorbing_states, ...
	'debug', debug, ...
	'algorithm', 'ttexpsumst');


% end

