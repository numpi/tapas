function acc_times = test_examplest(N, sd)
%TEST_EXAMPLES

rng(sd);

% Number of test to run
k = 1;

acc_times = [];

for j = 1 : k
    acc_times = [ acc_times, input_n_mttf(N, 'ttexpsumst', ...
				  full((eye(N)+sprand(N,N,.5/N)) > 0)) ];
end

fprintf('Average solution time: %f secs -- variance: %f \n', mean(acc_times), ...
    std(acc_times));

end

