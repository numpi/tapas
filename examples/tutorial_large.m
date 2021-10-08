function [MTTF, CMTTF, CMTTB] = tutorial_large(k)

%%
% setup of both TT-Toolbox and TAPAS
%%
cd ../../TT-Toolbox/
setup
cd ../tapas/
setup
cd examples/
%%
% define the infinitesimal generator matrix implicitly
%%
n = 2+k; % number of submodels
%%
% set parameters
l1 = 8.333e-2;
d2 = 8.333e-2;
m1 = 8.333e-2;%0.5;
m2 = 8.333e-2;%0.5;
m3 = 8.333e-2;%0.5;
d1 = 1;
t3 = 1;
%%%
sz = [3*ones(1,n-1) 4]; % sizes of the submodels
R = cell(1, n);
W = cell(2, n);
for j = 1 : k
    R{j} = tt_matrix([0   l1 0; ...
                      m1  0  0; ...
                      0   0  0]);
end
R{k+1} = tt_matrix([0   0  0; ...
                  m2  0  0;
                  0   0  0]);
R{k+2} = tt_matrix([0   0  0 0; ...
                  m3  0  0 t3; ...
                  0   0  0 0; ...
                  0   0  0 0]);

for j = 1 : k
    W{1,j} = tt_matrix([0  0  0; ...
                        0  0  d1; ...
                        0  0  0]);
end
W{1,k+1} = tt_matrix([0  0  0; ...
                    0  0  1; ...
                    0  0  0]);
W{1,k+2} = tt_matrix([0  0  1  0; ...
                    0  0  1  0; ...
                    0  0  0  0; ...
                    0  0  0  1]);

for j = 1 : k
    W{2,j} = tt_matrix([1  0  0; ...
                        0  1  0; ...
                        0  0  1]);
end
W{2,k+1} = tt_matrix([0  d2  0; ...
                    0  0   0; ...
                    0  0   0]);
W{2,k+2} = tt_matrix([0  1  0  0; ...
                    0  0  0  0; ...
                    0  0  0  0; ...
                    0  0  0  1]);
%%
% define absorbing states of interest
%%
F = [3*ones(1,n-1) 3];
B = [3*ones(1,n-1) 4];
absorbing_states =[F; B];
%%
% define initial probability and reward vectors
%%
pi0 = ktt_ej(sz, ones(1, n));
r   = round(ktt_ones(sz) - ktt_ej(sz, F) - ktt_ej(sz, B), 1.e-8);
%%
% set tolerances
%%
ttol = 1e-12;
tol = 1e-4;
%%
% evaluate the measures
%%

method = 'amen';%'spantree'; %'gmres';

MTTF = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
				    'algorithm', method, 'ttol', ttol,  ...
				    'absorbing_states', absorbing_states, 'tol', tol);                
[CMTTF, time] = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
		 			         'algorithm', method, ...
    				          'absorbing_states', absorbing_states, ...
					          'conditional_indices', F, ...
					          'ttol', ttol, 'tol', tol);
[CMTTB, time] = eval_measure('cond_etta', pi0, r, R, W, 'debug', true, ...
		 			         'algorithm', method, ...
    				          'absorbing_states', absorbing_states, ...
					          'conditional_indices', B, ...
					          'ttol', ttol, 'tol', tol);
end