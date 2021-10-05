function T = createTopology(n, p, graph)

switch graph

case 'star'

    T = eye(n) + eye(n, 1) * ones(1, n);

    for i = 1: n
      if rand < p
        % Create a second connection on this row
        j = randi(n);
        T(i,j) = T(i,j) + 1;
      end
    end
    
case 'starnoloops'
    T = createTopology(n, p, 'star');
    
    % Remove self loops of length 2
    for i = 1 : n
        for j = i+1 : n
            if T(i,j) && T(j,i)
                T(j,i) = 0;
            end
        end
    end

case 'erdos'
    T = eye(n) + ( rand(n,n) < 1/(n-1) );

    for i = 1 : n
        idx = find(T(:,i) > 0);
        if length(idx) == 1
            j = randi(n-1);
            T(j + (j >= i), i) = 1;
        end
        idx = find(T(i, :) > 0);
        if length(idx) > 3
            T(i,:) = 0;
            T(i,idx(1:2)) = 1;
            T(i,i) = 1;
        end
   end

case 'randperm'

    T = eye(n, n);

    T = T + T(randperm(n), :);

    for i = 1: n
      if rand < p
        % Create a second connection on this row
        j = randi(n);
        T(i,j) = T(i,j) + 1;
      end
	end
	

case 'rand'

    T = eye(n, n);

    for i = 1: n
      if rand < p
        % Create a second connection on this row
        j = randi(n);
        T(i,j) = T(i,j) + 1;
      end
    end	

case 'cycle'

    T = eye(n);
    T = T + T([2:n,1],:);

    for i = 1 : n
       if rand < p
          j = randi(n-2);
          T(i,j + 2*(j>=i)) = 1;
       end
    end

case 'tree'
    G = digraph;

    lvl = ceil(log2(n));

    for j = 2 : lvl
        counter = 0;

        for i = 2^(j-2) : 2^(j-1)-1
            % Connect each node with its two children
            G = addedge(G, i, 2^(j-1)+counter);
            counter = counter + 1;

            if counter + 2^(j-1) == n+1
                break;
            end

            G = addedge(G, i, 2^(j-1)+counter);
            counter = counter + 1;

            if counter + 2^(j-1) == n+1
                break;
            end
        end
    end
    T = full(adjacency(G)) + eye(n);

otherwise
    error('Unsupported graph');

end

T = (T > 0);

end
