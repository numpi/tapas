function F = ktt_format(A)
%KTT_FORMAT 

if isa(A, 'tt_matrix') || isa(A, 'tt_tensor')
	F = 'tt';
	return;
end

if issparse(A)
	F = 'sparse';
	return;
end

error('Unsupported format in A');

end

