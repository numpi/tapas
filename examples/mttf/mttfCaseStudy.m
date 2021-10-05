function [R, W, M] = mttfCaseStudy(n, topology, lambdah, lambdahp, lambdas, lambdasl, mu)
%LARGEMODEL produces matrices R and W and markings M
%   INPUTS: n:        number of system components,
%           topology: adjacency matrix of the attack propagation graph
%           lambdah:  vector of constant rates for the hardware failure
%           lambdahp: vector of constant rates for the hardware failure
%                     when the software has failed.
%           lambdas:  vector of constant rates for the software failure
%           lambdasl: vector of constant rates of local failure for the software 
%           mu:       recovery rates

% Check dimension
if size(topology) ~= [n, n]
    fprintf('ERROR: topology is required to be a %dx%d matrix\n',n);
    return;
end

R = cell(1, n);
W = cell(n, n);
M = cell(1,n);

for i = 1 : n
%     R{i} = zeros(3,3);
%     R{i}(1:2,3) = [ lambdah(i); lambdahp(i) ];
    R{i} = [0     lambdasl(i) lambdah(i);
            mu(i) 0           lambdahp(i);
            0     0           0];
    M{i} = [1 1 ;
            1 0 ;
            0 0];
end

for i = 1 : n
    for j = 1 : n
        if topology(i,j)
            % if there are more than 1 exterior attacks
            W{i,j} = zeros(3,3);
            
            if i == j % flag == 0
                W{i,j}(1,2) = lambdas(i);
            else
                W{i,j} = [ 0 1 0 ; 0 1 0 ; 0 0 1 ];
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
