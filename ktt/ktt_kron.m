function C = ktt_kron(varargin)
%

k = length(varargin);

if isa(varargin{1}, 'tt_matrix') || isa(varargin{1}, 'tt_tensor')
	C = varargin{1};
	for j = 2 : k
		C = tkron(varargin{j}, C);
	end
else
	C = varargin{1};
	for j = 2 : k
		C = kron(C, varargin{j});
	end
end

end

