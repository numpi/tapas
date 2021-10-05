function [m, time, y] = inv_tt_regular_splitting(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_TT_REGULAR_SPLITTING 

timer = tic;
        
% DA are the factors of the Kronecker sum for (gamma * eye(n) - R)
% Note that the sign is changed with respect to R - gamma*eye. 
DA = inv_computeDA(R, W, shift);
[~, Delta] = ktt_infgen(R, W, ttol);
[~, ~, Deltap, Wsync] = inv_computeS(R, W, Delta, absorbing_states, shift, ttol);

n = arrayfun(@(i) size(R{i}, 1), 1 : length(R));

% FIXME: we need to evaluate the number of terms in a sensible way.
expn = ceil(- 6 * log(ttol)/pi);
expsums_method = 'sinc';
maxsteps = inf;
debug = true;

en = ktt_ej(n, n, 'tt');

% In case we need to debug
%M = full(ktt_kronsum(R{:})) + full(Deltap); M(:,end) = M(:,end) - full(ktt_kronsum(R{:})) * full(en);
%S = zeros(prod(n)); S(:,end) = full(ktt_kronsum(R{:}) + Wsync + Delta) * full(en); S(end,end) = S(end,end) + 1;
%N = -full(Delta) + full(Deltap) - full(Wsync) + S; N(:,end) = N(:,end) - full(ktt_kronsum(R{:})) * full(en);

% Precompute f                
% f = (R - gamma * eye(prod(n))) \ R(:,end);
Rend = round(ktt_kronsum(R{:}) * en, ttol);
f = ttexpsummldivide(DA, -Rend, expn, ttol, expsums_method);        
g = f * (1./ (1 - dot(en, f)));

Su = round(Rend + (Wsync*en) + (Delta*en), ttol);

Mb = ttexpsummldivide(DA, -r, expn, ttol, expsums_method); 
Mb = round(Mb - g * dot(en, Mb), ttol);

x0 = Mb;
x = x0;

nrmx0 = norm(Mb);

j = 0;

while j < maxsteps
    j = j + 1;

    % Adjust the truncation tolerance based on the accuracy that we
    % have achieved as of now
    ltol = ttol; % max(ttol, rel_change / cnd / 10);

    xold = x;
    % x = M \ (N * x) + x0;

    % Compute l = N * x
    l = -Wsync * x;
    l = round(l - Delta * x + Deltap * x, ltol);
    enx = dot(en, x);
    l = round(l + Su * enx - Rend * enx, ltol);

    % Solve the linear system M*x = l
    x = ttexpsummldivide(DA, -l, expn, ltol, expsums_method); 
    x = round(x - g * dot(en, x), ltol);

    % Estimate the spectral radius            
    rho = norm(x) / norm(xold);

    % Update the iterate
    x = x + Mb;


    res_a = inf;

    % Estimate for the error
    if rho < 1
        err_est = nrmx0 / norm(x) * rho^(j+1) / (1 - rho);
    else
        err_est = inf;
    end            

    if debug
        fprintf('Step %d, err. est. = %e, m = %e, ranks = %d, rho = %e\n', ...
                j, err_est, -dot(pi0, x), max(rank(x)), rho);
    end

    if err_est < tol
        break;
    end

    if res_a < tol
        x = xa;
        break;
    end
end

m = -dot(pi0, x);

y = x;

t = toc(timer);
time = t;

end

