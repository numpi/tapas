function M = ktt_outerprod(V, W)
%KTT_OUTERPROD 

fmt = ktt_format(V);

m = size(V);
n = size(W);

if isa(V, 'tt_matrix') || isa(W, 'tt_matrix')
	error('tt_matrix format unsupported in ktt_outerprod');
end

switch fmt
	case 'tt'
		M = kron(W, V); M = tt_matrix(M, m, n);
	case 'sparse'
		M = V * W';
	otherwise
		error('Unsupported format');
end

end

