%%
% setup of both TT-Toolbox and TAPAS
%%
cd ../../TT-Toolbox/
setup
cd ../tapas/
setup
cd examples/
%%
% set parameters
%%
l1 = 0.5;
l2 = 0.5;
m1 = 1;
m2 = 1;
m3 = 1;
d1 = 0.01;
t3 = 0.1;
%%
% define the infinitesimal generator matrix implicitly
%%
n = 3; % number of submodels
sz = [3 3 4]; % sizes of the submodels
R = cell(1, n);
W = cell(2, n);
R{1} = tt_matrix([0   l1 0; ...
                  m1  0  0; ...
                  0   0  0]);
R{2} = tt_matrix([0   0  0; ...
                  m2  0  0;
                  0   0  0]);
R{3} = tt_matrix([0   0  0 0; ...
                  m3  0  0 t3; ...
                  0   0  0 0; ...
                  0   0  0 0]);

W{1,1} = tt_matrix([0  0  0; ...
                    0  0  d1; ...
                    0  0  0]);
W{1,2} = tt_matrix([0  0  0; ...
                    0  0  1; ...
                    0  0  0]);
W{1,3} = tt_matrix([0  0  1  0; ...
                    0  0  1  0; ...
                    0  0  0  0; ...
                    0  0  0  1]);

W{2,1} = tt_matrix([1  0  0; ...
                    0  1  0; ...
                    0  0  1]);
W{2,2} = tt_matrix([0  l2  0; ...
                    0  0   0; ...
                    0  0   0]);
W{2,3} = tt_matrix([0  1  0  0; ...
                    0  0  0  0; ...
                    0  0  0  0; ...
                    0  0  0  1]);
%%
% define absorbing states of interest
%%
F = [3 3 3];
B = [3 3 4];
absorbing_states =[F; B];
%%
% define initial probability and reward vectors
%%
pi0 = ktt_ej(sz, ones(1, n));
r   = round(ktt_ones(sz) - ktt_ej(sz, F1) - ktt_ej(sz, F2), 1.e-8);
%%
% evaluate the measures
%%
m = eval_measure('inv', pi0, r, R, W, 'debug', true, ...
				 'algorithm', 'spantree', 'ttol', 1e-12,  ...
				 'absorbing_states', absorbing_states);