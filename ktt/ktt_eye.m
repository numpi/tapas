function I = ktt_eye(n, fmt)
%KTT_EYE 

if ~exist('fmt', 'var')
	fmt = 'tt';
end

if isempty(n)
	error('Empty matrices are not (yet) supported');
end

switch fmt
	case 'tt'
		I = tt_eye(n(end:-1:1));
	case 'sparse'
		I = speye(prod(n));
	otherwise
		error('Unsupported format');
end

end

