function [S, DeltapC, Deltap, Wsync] = inv_computeS(R, W, Delta, absorbing_states, shift, ttol)
%INV_COMPUTES 

fmt = ktt_format(R{1});
n = arrayfun(@(i) size(R{i}, 1), 1 : length(R));

DeltapC = diagblocks(R, W, shift);
Deltap = -round(ktt_kronsum(DeltapC{:}), ttol);

Wsync = ktt_zeros(n, n, fmt);
for i = 1 : size(W, 1)
	Wsync = round(Wsync + ktt_kron(W{i,:}), ttol);
end

A1 = round(ktt_kronsum(R{:}), ttol);
A2 = round(Wsync + (Delta - Deltap), ttol);

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

end

