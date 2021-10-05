function time = example_BCfailure_condetta(method)
%
% Constructs the example decribed in the paper "New Method for 
% the Evaluation of (Conditional) Time To Failure and its Moments
% Trough Implicit Reward Structure", by G. Masetti, L. Robol, S. Chiaradonna,
%   F. Di Giandomenico, submitted.

nchanges = 0;% number of different value of pC between 0 and 1-pEP

n = 6;% number of components

% {{true,false,false,false,false,false},
% {true,true,false,false,false,false},
% {true,false,false,false,false,false},
% {true,true,false,false,false,false},
% {true,false,false,false,false,false},
% {true,false,false,false,false,false}}
% isDelta = [1 0 0 0 0 0;
%             1 1 0 0 0 0;
%             1 0 1 0 0 0;
%             1 1 0 1 0 0;
%             1 0 0 0 1 0;
%             1 0 0 0 0 1];
topology = [1 1 1 1 1 1;
            0 0 0 1 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0];
        
if ~exist('topology', 'var')
    topology = createTopology(n, 0.2, 'starnoloops');    
end

pp = symrcm(topology);
disp(topology(pp,pp));

% n = 2;
% topology = [1 1
%             0 0];

lambdaB = 18*ones(n,1)

lambdaD = 10*ones(n,1)

lambdaW = 342*ones(n,1)

lambdaE = ones(n,1)

lambdaEP = 5*ones(n,1)

lambdaC = 0.041665*ones(n,1)
  
% Construct Rs and Ws
[R, W] = BCfailure(n, topology(pp,pp), lambdaB(pp), lambdaC(pp), lambdaD(pp), ...
	lambdaW(pp), lambdaE(pp), lambdaEP(pp));

P = eye(5); 
P(1:3,1:3) = 0;
w = cell(1, n);
for j = 1 : n
	w{j} = tt_matrix(P);
end
absorbing_states = ktt_kron(w{:});

pi0 = ktt_ej(5*ones(1,n), ones(1,n));

v = [ 0 ; 0 ; 0 ; 1 ; 1 ];
w = cell(1, n); for j = 1 : n; w{j} = tt_tensor(v); end

r = round(ktt_ones(5*ones(1,n)) - ktt_kron(w{:}), 1e-8);

% C is 4-th, B is 5-th, state
% conditional_indices = 4 * ones(1, n);
% conditional_indices = 5 * ones(1, n);
conditional_indices = [ 5 * ones(1, n); 4 * ones(1, n) ];

% conditional_indices = [4 5];
% conditional_indices = [5 4];

% % Compute the measure
% if n <= 6
% 	m = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
% 				 'algorithm', 'spantree', ...
% 				 'conditional_indices', conditional_indices, ...
% 				 'absorbing_states', absorbing_states);
% end

shift = 0;
it_mult = false;

if strcmp(method, 'ttexpsums2')
	shift = 2e6;
end

[m, time] = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
					   'algorithm', method, 'shift', shift, ...
					   'absorbing_states', absorbing_states, ...
					   'conditional_indices', conditional_indices, ...
					   'ttol', 1e-10, 'tol', 1e-4, ...
					   'iterative_mult', it_mult, 'use_sinc', false, ...
					   'interval_report', 10);

end

