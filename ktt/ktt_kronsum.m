function C = ktt_kronsum(varargin)
%KTT_KRONSUM 

k = length(varargin);
n = zeros(1, k);

fmt = ktt_format(varargin{1});

for j = 1 : k
	n(j) = size(varargin{j}, 1);
end

if k == 1
	C = varargin{1};
	return;
end

C = ktt_kron(varargin{1}, ktt_eye(n(2:end), fmt)) + ...
	ktt_kron(ktt_eye(n(1:end-1), fmt), varargin{end});

for j = 2 : k-1
	C = C + ktt_kron(ktt_eye(n(1:j-1), fmt), ...
					 varargin{j}, ...
					 ktt_eye(n(j+1:end), fmt));
end

end

