function r = lanczos_fa(f, pi0, R, W, r, absorbing_states)
%LANCZOS_FA Evaluate w' * f(A) * v

% Phase 1: prepare the matrix Q that defines the matrix function
Wsync = ktt_zeros(n, n, fmt);
for i = 1 : size(W, 1)
	Wsync = round(Wsync + ktt_kron(W{i,:}), ttol);
end

QQ = round(ktt_kronsum(R{:}) + Wsync, ttol);

A1 = round(ktt_kronsum(R{:}), ttol);

% Construct the low-rank correction to implicitly deflate the absorbing states
if isa(absorbing_states, 'tt_matrix')
	S = absorbing_states;
	S = round(A1 * S, ttol) + round(A2 * S, ttol);
	S = round(S, ttol);
else
	S = ktt_zeros(n, n, fmt);
	for i = 1 : size(absorbing_states, 1)
		ej = ktt_ej(n, absorbing_states(i,:), fmt);
		S = round(S + ktt_outerprod(A1*ej + A2*ej, ej), ttol);
	end
end

Q = round(QQ + Delta - S, ttol);

Q = full(Q);

% Prepare the vectors v, w
V = [pi0, r];
W = [r, pi0];



end

