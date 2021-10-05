function [R, W, M] = BCfailure(n, topology, lambdaB, lambdaC, lambdaD, lambdaW, lambdaE, lambdaEP, l)
%BCFAILURE produces matrices R and W and markings M
%   INPUTS: n:        number of system components,
%           topology: adjacency matrix of the failure propagation graph
%           lambdaB:  benign component failure rates
%           lambdaC:  catastrofic component failure rate
%           lambdaD:  component failure detection rate
%           lambdaW:  component recovery rate
%           lambdaEP: failure propagation rate
%
% The model constructed is described in the paper "New Method for 
% the Evaluation of (Conditional) Time To Failure and its Moments
% Trough Implicit Reward Structure", by G. Masetti, L. Robol, S. Chiaradonna,
% F. Di Giandomenico, submitted.

% Check dimension
if size(topology) ~= [n, n]
    fprintf('ERROR: topology is required to be a %dx%d matrix\n',n);
    return;
end

R = cell(1, n);
W = cell(n, n);
M = cell(1,n);

for i = 1 : n
    if i <= l
        R{i} = zeros(4,4);
        R{i}(1,2) = lambdaE(i);
        R{i}(2,3) = lambdaD(i);
        R{i}(3,1) = lambdaW(i);
        R{i}(2,4) = lambdaC(i);
        M{i} = eye(4,4);
    else        
        R{i} = zeros(5,5);
        R{i}(1,2) = lambdaE(i);
        R{i}(2,3) = lambdaD(i);
        R{i}(3,1) = lambdaW(i);
        R{i}(3,5) = lambdaB(i);
        R{i}(2,4) = lambdaC(i);
        M{i} = eye(5,5);
    end
end

% manage TEP
for i = 1 : n
    for j = 1 : n
        if topology(i,j)
            if j <= l
                W{i,j} = zeros(4,4);
            else
                W{i,j} = zeros(5,5);
            end
            
            if i == j
                W{i,j}(2,2) = lambdaEP(i);
            else
                if j <= l
                    W{i,j} = eye(4,4);
                else
                    W{i,j} = eye(5,5);
                end
				W{i,j}(1:3,1:3) = 0; 
				W{i,j}(1:3,2) = 1;
            end
        else
            % if component j is subject to NO exterior attack
            if j <= l
                W{i,j} = eye(4,4);
            else
                W{i,j} = eye(5,5);
            end
        end
    end
end

% manage TC
%for i = n+1 : 2*n
%    for j = 1 : n
%            W{i,j} = zeros(5,5);
%           if i-n == j
%               W{i,j}(2,4) = lambdaC(i-n);
%            else
%                W{i,j} = eye(5,5);
%				W{i,j}(1:3,1:3) = 0; 
%				W{i,j}(1:3,4) = 1;
%            end
%    end
%end

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
