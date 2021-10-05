function acc_times = test_examples_itmult(N, sd)
%TEST_EXAMPLES

rng(sd);

% Number of test to run
k = 1;

acc_times = [];

for j = 1 : k
    acc_times = [ acc_times, input_n_density01(N, 'ttexpsums2', createTopology(N, 0.2, 'starnoloops'), true) ];
end

fprintf('Average solution time: %f secs -- variance: %f \n', mean(acc_times), ...
    std(acc_times));

end

