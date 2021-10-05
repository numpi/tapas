function [R, W, M] = CCFModel2(n, topology, iso, cas, mu, ok, red)
%CCFMODEL Common-Cause-of-Failure
%
% This model describes a system of N components, each can fail
% independently with rate ISO(j), j = 1, ..., N; alternatively, the
% component J can cause a failure to itself and all the components connected 
% according the given TOPOLOGY. 
%
% Each failed component is repaired with rate MU(j). 

% Count the rows of T where we need to have some synchronization
v = ( topology - eye(n) ) * ones(n, 1); v = find(v ~= 0);

R = cell(1, n);
W = cell(length(v), n);
M = cell(1,n);

% States are as follows:
%
% i \in { 1 , ..., red } : red-i+1 replicas are working
% i == red + 1           : Ok
% i == red + 2           : Failure
%
% Transitions:
%   i -> i + 1 for i { 1, ..., red-1 } with rates iso(J) (local)
%   red -> red + 2 with rate iso(J) (local)
%   { 1, ..., red } => red + 2 rate cas(J) (sync) 
%   i -> i - 1 for i { 2, ..., red } with rates mu(J) (local)
%   red+2 -> red (local*)
%                *ONLY ENABLED UNLESS ALL ARE FAILED!
%
%   { 1, ..., red } -> red + 1 (Ok, local) with rate ok(J)

for j = 1 : n
    R{j} = zeros(red+2,red+2);
	
	for i = 1 : red-1
		R{j}(i,i+1) = (red-i+1) * iso(j);
	end
	
	R{j}(red, red+2) = iso(j);
	
	for i = 2 : red
		R{j}(i,i-1) = mu(j);
	end
	
	R{j}(red+2,red) = mu(j);
	
	for i = 1 : red
		R{j}(i,red+1) = ok(j);
	end
end

for i = 1 : length(v)
   for j = 1 : n
      if ( topology(v(i),j) == 1 )
          W{i,j} = zeros(red+2,red+2);
          W{i,j}(1:red,red+2) = cas(v(i));
      else
          % no interaction between i and j
    	  W{i,j} = eye(red+2,red+2);
      end
   end
end


for j = 1 : n
	W{end+1,j} = zeros(red+2,red+2); W{end,j}(red+2,red) = -mu(j);
	for i = 1 : n
		if i ~= j
			W{end,i} = zeros(red+2,red+2); W{end,i}(red+2,red+2) = 1;
		end
	end
end

end
