function acc_times = test_examples(N, sd)
%TEST_EXAMPLES

rng(sd);

% Number of test to run
k = 1;

acc_times = [];

for j = 1 : k
    acc_times = [ acc_times, input_n_density01(N, 'ttexpsums2') ];
end

fprintf('Average solution time: %f secs -- variance: %f \n', mean(acc_times), ...
    std(acc_times));

end

