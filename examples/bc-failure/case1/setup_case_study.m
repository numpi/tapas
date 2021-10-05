% Topology
rng(2)

% Tolerances
ttol = 1e-10;
tol = 1e-4;

%  topology = full((eye(n,n)+sprand(n,n, 1 / n)) > 0)
if ~exist('topology', 'var')
    topology = createTopology(n, 0.2, 'starnoloops');
end

% Topology generate with Erdos plus a superdiagonal correction equal to all
% ones. 
load('topology', 'T');
topology = T(1:n, 1:n);
xx = rand(n-1,1) > .5;
% yy = rand(n-1,1) > .5;
topology = eye(n) + diag(xx .* ones(n-1,1),1) + diag(ones(n-1,1),-1);
topology = topology > 0;

disp(topology);

% Number of systems with only one absorbing state
l = ceil(n/2);

% Construct Rs and Ws
[R, W] = BCfailure(n, topology, lambdaB, lambdaC, lambdaD, ...
	lambdaW, lambdaE, lambdaEP, l);

P = eye(5); P(1:3,1:3) = 0;
w = cell(1, n);
for j = 1 : n
    if j <= l
        w{j} = tt_matrix(P(1:4,1:4));
    else
        w{j} = tt_matrix(P);
    end
end
absorbing_states = ktt_kron(w{:});

pi0 = ktt_ej([ 4*ones(1,l), 5*ones(1,n-l) ], ones(1,n));

v = [ 0 ; 0 ; 0 ; 1 ; 1 ];
w = cell(1, n); for j = 1 : n
    if j <= l
        w{j} = tt_tensor(v(1:4));
    else
        w{j} = tt_tensor(v); 
    end
end

r = round(ktt_ones([ 4*ones(1,l), 5*ones(1,n-l) ]) - ktt_kron(w{:}), 1e-8);
