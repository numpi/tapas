function [R, W, M] = largeModel(n, topology, lambda, mu, p)
%LARGEMODEL produces matrices R and W and markings M
%   INPUTS: n:        number of system components,
%           topology: adjacency matrix of the attack propagation graph
%           lambda:   vector of constant rates for the attack
%           mu:       vector of constant rates for the recovery
%           p:        vector of constant probabilities
%
% The model constructed is described in the paper "Stochastic modeling and 
% evaluation of large interdependent composed models through Kronecker 
% algebra and Exponential sums", by G. Masetti, L. Robol, S. Chiaradonna,
% F. Di Giandomenico, submitted.

% Check dimension
if size(topology) ~= [n, n]
    fprintf('ERROR: topology is required to be a %dx%d matrix\n',n);
    return;
end

R = cell(1, n);
W = cell(n, n);
M = cell(1,n);

% Stochastic Petri Net of component i:
%
% Places:
% Safe_i, Tacken_i, CatastroficTaken_i.
% Local transitions:
% LocalAttack_i, LocalRecovery_i.
% Synchronization transitions:
% SyncAttack_j.
%
% Arcs:
% Safe_i -> LocalAttack_i
% LocalAttack_i -> Taken_i
% Taken_i -> LocalRecovery_i
% LocalRecovery_i -> Safe_i
% Safe_i -> SyncAttack_j

for i = 1 : n
    R{i} = zeros(3,3);
    R{i}(1,2) = p(i)*lambda(i);
    R{i}(2,1) = mu(i);
    M{i} = [1 0 0;
            0 1 0;
            0 0 1];
end

for i = 1 : n
    % numberOfDeps = sum(topology(:,i));
    flag = 0;
    for j = 1 : n
        if topology(i,j)
            % if there are more than 1 exterior attacks
            W{i,j} = zeros(3,3);
            
            if i == j % flag == 0
                W{i,j}(1,3) = lambda(i)*(1-p(i));
            else
                W{i,j}(1:3,3) = 1;                
            end
        else
            % if component i is subject to NO exterior attack
            W{i,j} = eye(3,3);
        end
    end
end

% convert Rs and Ws in TT-format
for i = 1 : n    
    R{i} = tt_matrix(R{i});
end

for i = 1 : n
    for j = 1 : n
        W{i,j} = tt_matrix(W{i,j});
    end
end

end
