function time = input_n_density01_BCfailure_condetta(n, method, casenumber)
%
% Constructs the example decribed in the paper "New Method for 
% the Evaluation of (Conditional) Time To Failure and its Moments
% Trough Implicit Reward Structure", by G. Masetti, L. Robol, S. Chiaradonna,
%   F. Di Giandomenico, submitted.

lambdaB_min = 18;
lambdaB_max = 20;
lambdaB = lambdaB_min + (lambdaB_max-lambdaB_min)*rand(n,1)

lambdaC_min = 0.041665-0.01;
lambdaC_max = 0.041665+0.01;
lambdaC = lambdaC_min + (lambdaC_max-lambdaC_min)*rand(n,1)

lambdaD_min = 9;
lambdaD_max = 11;
lambdaD = lambdaD_min + (lambdaD_max-lambdaD_min)*rand(n,1)

lambdaW_min = 320;
lambdaW_max = 350;
lambdaW = lambdaW_min + (lambdaW_max-lambdaW_min)*rand(n,1)

lambdaE_min = 1;
lambdaE_max = 1.5;
lambdaE = lambdaE_min + (lambdaE_max-lambdaE_min)*rand(n,1)

lambdaEP_min = 0.9;
lambdaEP_max = 1.1;
lambdaEP = lambdaEP_min + (lambdaEP_max-lambdaEP_min)*rand(n,1)

setup_case_study;

switch casenumber
    case 1
        % Catastrophic case
        conditional_indices = 4 * ones(1, n);
    case 2
        % All in the benign case
        conditional_indices = [ 4 * ones(1,l), 5 * ones(1, n-l) ];
    case 3
        % Catastrophic + exactly one benign case, for i > l
        conditional_indices = zeros(0, n);
        for i = l+1 : n
            w = 4 * ones(1, n); 
            w(i) = 5;
            conditional_indices = [ conditional_indices ; w ];
        end
    case 4
        % Catastrophic + exactly two benign cases, also for i > l
        conditional_indices = zeros(0, n);
        for i1 = l+1 : n
            for i2 = i1 + 1 : n
                w = 4 * ones(1, n); 
                w(i1) = 5;
                w(i2) = 5;
                conditional_indices = [ conditional_indices ; w ];
            end
        end
end

[m, time] = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
					   'algorithm', method, 'batch_size', 2, ...
					   'absorbing_states', absorbing_states, ...
					   'conditional_indices', conditional_indices, ...
					   'ttol', ttol, 'tol', tol);

end

