function Z = ktt_zeros(varargin)
%KTT_ZEROS 

ninputs = length(varargin);

if ischar(varargin{end})
	fmt = varargin{end};
	ninputs = ninputs - 1;
else
	fmt = 'tt';
end

if ninputs > 2
	error('Unsupported number of inputs');
end

switch fmt
	case 'tt'
		if ninputs == 1
			Z = tt_zeros(varargin{1}(end:-1:1));
		else
			Z = tt_matrix(tt_zeros(varargin{1}(end:-1:1) .* varargin{2}(end:-1:1)), varargin{1}(end:-1:1), varargin{2}(end:-1:1));
		end
	case 'sparse'
		m = varargin{1};
		if ninputs == 1
			Z = sparse(prod(m), 1);
		else
			n = varargin{2};
			Z = sparse(prod(m), prod(n));
		end
	otherwise
		error('Unsupported format');
end

end

