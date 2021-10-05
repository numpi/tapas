function time = example_varying_pC_BCfailure_mean_tta_and_variance(method)

nchanges = 0;% number of different value of pC between 0 and 1-pEP

n = 6;% number of components
% 
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

% n = 2;
% topology = [1 1
%             0 0];

disp(topology);

lambdaB = 18*ones(n,1)

lambdaD = 10*ones(n,1)

lambdaW = 342*ones(n,1)

lambdaE = ones(n,1)

lambdaEP = 5*ones(n,1)

pp = symrcm(topology);
disp(topology(pp,pp));

P = eye(5); P(1:3,1:3) = 0;
w = cell(1, n);
for j = 1 : n
	w{j} = tt_matrix(P);
end
absorbing_states = ktt_kron(w{:});

pi0 = ktt_ej(5*ones(1,n), ones(1,n));

v = [ 0 ; 0 ; 0 ; 1 ; 1 ];
w = cell(1, n); for j = 1 : n; w{j} = tt_tensor(v); end

r = round(ktt_ones(5*ones(1,n)) - ktt_kron(w{:}), 1e-8);

shift = 0.01;
it_mult = true;

if strcmp(method, 'ttexpsums2')
	shift = 1e6;
end

results = zeros(nchanges+1,3);

for k = 0 : nchanges
%     lambdaC_min = 0.1+k*0.05;
%     lambdaC_max = 0.2+k*0.05;
%     lambdaC = lambdaC_min + (lambdaC_max-lambdaC_min)*ones(n,1)
    lambdaC = 0.041665*ones(n,1)
    
    % Construct Rs and Ws
    [R, W] = BCfailure(n, topology(pp, pp), lambdaB(pp), lambdaC(pp), lambdaD(pp), ...
        lambdaW(pp), lambdaE(pp), lambdaEP(pp));
    
    % Debug
%     Q = infgen(R, W, 1e-8, true);
%     tt = find(sum(full(abs(Q)), 2) == 0);
%     tau = setdiff(1:size(Q,1), tt);
%     Qt = Q(tau, tau);
%     fpi0 = eye(size(Qt, 1), 1);
%     fr = ones(size(Qt, 1), 1);
%     m1 = - fpi0' * (Qt \ fr);
    
    [mtta, time] = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
                           'algorithm', method, 'shift', shift, ...
                           'absorbing_states', absorbing_states, ...
                           'ttol', 1e-10, 'tol', 1e-4, ...
                           'iterative_mult', it_mult, 'use_sinc', true, ...
                           'interval_report', 10);
    [var, time] = eval_measure('tta_variance', pi0, r, R, W, 'debug', true, ...
                           'algorithm', method, 'shift', shift, ...
                           'absorbing_states', absorbing_states, ...
                           'ttol', 1e-10, 'tol', 1e-4, ...
                           'iterative_mult', it_mult, 'use_sinc', true, ...
                           'interval_report', 10);
    results(k+1,:) = [lambdaC(1), mtta, var];
end

fprintf("lambdaC\t\tMTTF\t\tVar\n");
fprintf("%f\t%f\t%f\n", results');
                   
end

