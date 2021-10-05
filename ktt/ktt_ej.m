function ej = ktt_ej(n, j, fmt)
%KTT_EJ 

if ~exist('fmt', 'var')
	fmt = 'tt';
end

k = length(n);
v = cell(1, k);

j = j(:);
n = n(:);

if any(j > n) || any(j <= 0)
	error('Indices out of bounds');
end

switch fmt
	case 'tt'
		for i = 1 : k
			v{i} = zeros(n(i), 1); v{i}(j(i)) = 1;
			v{i} = tt_tensor(v{i});
		end
	case 'sparse'
		for i = 1 : k
			v{i} = spzeros(n(i), 1); v{i}(j(i)) = 1;
		end
	otherwise
		error('Unsupported format');
end

ej = ktt_kron(v{:});

end

