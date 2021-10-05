function [Q, Delta] = ktt_infgen(R, W, ttol)
%KTT_INFGEN 

k = length(R);
n = arrayfun(@(i) size(R{i}, 1), 1 : k);

fmt = ktt_format(R{1});

Q = ktt_kronsum(R{:}); 

for i = 1 : size(W, 1)
	Q = Q + ktt_kron(W{i,:});
    if strcmp(fmt, 'tt')
        Q = round(Q, ttol);
    end
end

d = Q * ktt_ones(n, fmt);

switch fmt
	case 'sparse'
        Delta = -spdiags(d, 0, prod(n), prod(n));
		Q = Q + Delta;
	case 'tt'
        Delta = -diag(d);
		Q = Q + Delta;
        Q = round(Q, ttol);
end

end

