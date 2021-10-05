function M = ktt_iterative_mult(X, ttol, debug)
%KTT_ITERATIVE_MULT Iteratively approximate X * X
%
% The method decomposes X = X0 + dX, and approximates X*X by
%
%  X * X = X0*X0 + X0*dX + dX*X0 + KTT_ITERATIVE_MULT(dX, dX, ...),
%
% where the recursive call is done with a larger threshold induced by the
% fact that dX is (hopefully) smaller than X.
%
% EXPERIMENTAL CODE! Might not converge.

if ~exist('debug', 'var')
	debug = false;
end

nrmX = norm(X);

ltol = ttol;

X0 = round(X, ttol, ceil(sqrt(length(size(X)))));
dX = round(X - X0, ttol);

M = round(X0 * X0, ttol);
M = round(M + dX * X0, ttol);
M = round(M + X0 * dX, ttol);

X = dX;

maxit = 10000;

maxrank = 8; % length(size(X)); % max(5, ceil(sqrt(length(size(X)))));

for j = 1 : maxit
	nrmdX = norm(dX);
	
	% Adjust the relative tolerance --we never allow to go below 1e-1
	rtol = min(1, ttol * (nrmX / nrmdX)^2);
	
	if debug
		fprintf('KTT_ITERATIVE_MULT :: rank(dX) = %d, rank(M) = %d, res = %e\n', ...
			max(rank(dX)), max(rank(M)), nrmdX^2 / nrmX^2);
	end	
	
	if nrmdX^2 < nrmX^2 * ttol * 100
		break;
	end
	
	X0 = round(X, rtol, maxrank);
	dX = round(X - X0, rtol);

	M = round(M + X0 * X0, ltol);
	M = round(M + dX * X0, ltol);
	M = round(M + X0 * dX, ltol);	
	
	X = dX;
end

end

