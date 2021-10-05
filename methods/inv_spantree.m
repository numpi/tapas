function [m, time, y] = inv_spantree(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_SPANTREE 

tspantree = tic;

debug = true;

k = length(R);
n = arrayfun(@(i) size(R{i}, 1), 1 : k);

% Construction of the matrix Q in sparse format
SR = cell(1, k); SW = cell(size(W,1), k);
for i = 1 : k
    SR{i} = sparse(full(R{i}));
    for j = 1 : size(W, 1)
        SW{j,i} = sparse(full(W{j,i}));
    end
end
Q = ktt_kronsum(SR{:}); 
for i = 1 : size(SW, 1)
    Q = Q + ktt_kron(SW{i,:});
end
d = Q * ktt_ones(n, 'sparse'); Q = Q - spdiags(d, 0, prod(n), prod(n));

% Find the reachable states by walking on the graph
G = digraph(abs(Q) > 1e-5, 'omitselfloops');
G = shortestpathtree(G, 1);
G = graph(G.Edges, G.Nodes);
t = minspantree(G, 'Root', 1);
idx = unique(t.Edges.EndNodes);

if debug
    fprintf(' - Number of reachable states: %d\n', length(idx));
end

% Construct dense representation of the rewards and of the pi0
fr = full(r); fpi0 = full(pi0);

% Construct the set of indices of the absorbing states
if isa(absorbing_states, 'tt_matrix')
    S = full(absorbing_states);
    abs_idx = find(sum(S, 2) > 1e-8);
else
    abs_idx = [];
    for j = 1 : size(absorbing_states, 1)
        nn = cumprod(n, 'reverse');
        abs_idx = [abs_idx, absorbing_states(j,end) + ...
                   sum((absorbing_states(j,1:end-1)-1).*nn(2:end))];
    end
end

idx = setdiff(idx, abs_idx);
y = Q(idx,idx) \ fr(idx); m = -dot(y, fpi0(idx));

% Convert to tt_tensor for compatibility
yy = zeros(length(fr), 1); yy(idx) = y;
y = tt_tensor(yy, size(r));

t = toc(tspantree);
time = t;

end

