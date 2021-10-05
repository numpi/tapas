function time = lambdaC_BCfailure_mean_tta_and_variance(method, k)

n = 6;% number of components

lambdaB_min = 18;
lambdaB_max = 20;
lambdaB = lambdaB_min + (lambdaB_max-lambdaB_min)*rand(n,1)

lambdaD_min = 9;
lambdaD_max = 11;
lambdaD = lambdaD_min + (lambdaD_max-lambdaD_min)*rand(n,1)

% lambdaW=342;
lambdaW_min = 320;
lambdaW_max = 350;
lambdaW = lambdaW_min + (lambdaW_max-lambdaW_min)*rand(n,1)

lambdaE_min = 1;
lambdaE_max = 1.5;
lambdaE = lambdaE_min + (lambdaE_max-lambdaE_min)*rand(n,1)

lambdaEP_min = 0.9;
lambdaEP_max = 1.1;
lambdaEP = lambdaEP_min + (lambdaEP_max-lambdaEP_min)*rand(n,1)

lambdaC = (k+1)*0.041665*ones(n,1);
    
setup_case_study;
        
[mtta, timemtta] = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol);
[var, timevar] = eval_measure('tta_variance', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, ...
                       'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol);
                   
results = [lambdaC(1), mtta, var, timemtta, timevar];

fprintf("lambdaC\t\tMTTA\t\tVar\t\tTimeMTTA\tTimeVAR\n");
fprintf("%f\t%f\t%f\t%f\t%f\n", results');

time = timemtta + timevar;
                   
end

