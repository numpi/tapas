function time = input_n_density01_BCfailure_moments_tta(n, method, topology)
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

[m1, timem1] = eval_measure('momentk', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, 'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol, 'moment', 1);
[m2, timem2] = eval_measure('momentk', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, 'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol, 'moment', 2);
[m3, timem3] = eval_measure('momentk', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, 'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol, 'moment', 3);
[m4, timem4] = eval_measure('momentk', pi0, r, R, W, 'debug', true, ...
                       'algorithm', method, 'absorbing_states', absorbing_states, ...
                       'ttol', ttol, 'tol', tol, 'moment', 4);
        

results = [lambdaC(1), m1, m2, m3, m4, timem1, timem2, timem3, timem4];

fprintf("lambdaC, m1, m2, m3, m4, TimeM1, TimeM2, TimeM3, TimeM4\n");
fprintf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", results');

skw = (m3-3*m2*m1-2*m1^3)/(m2-m1^2)^(3.0/2);

fprintf("skewness is %f\n", skw);

kur = (m4-4*m3*m1+6*m2*m1^2-3*m1^4)/(m2^2-2*m1*m2+m1^2);

fprintf("kurtosis is %f\n", kur);

time = timem4;

end

