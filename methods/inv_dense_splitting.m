function [m, time, y] = inv_dense_splitting(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_DENSE_SPLITTING 

timer = tic;

n = arrayfun(@(i) size(R{i}, 1), 1 : length(R));
maxsteps = inf;
debug = true;

[~, Delta] = ktt_infgen(R, W, ttol);
[~, ~, Deltap, Wsync] = inv_computeS(R, W, Delta, absorbing_states, shift, ttol);

% Q = full(infgen(R,W,ttol,true));
Delta = full(Delta);
R = full(ktt_kronsum(R{:}));
W = full(Wsync);        
Deltap = full(Deltap);        

% gamma = max(1, norm(Delta, inf)) + 1;
% Deltap = -gamma * eye(prod(n));

S = zeros(prod(n)); S(:,end) = R(:,end) + W(:,end) + Delta(:,end);
S(end,end) = S(end,end) + 1;

% Start the iteration
% M = R - gamma * eye(prod(n)); M(:,end) = M(:,end) - R(:,end);
M = R + Deltap; M(:,end) = M(:,end) - R(:,end);
% N = -Delta - gamma * eye(prod(n)) - W + S; N(:,end) = N(:,end) - R(:,end);
N = -Delta + Deltap - W + S; N(:,end) = N(:,end) - R(:,end);

fprintf('spectral radius=1-%e\n',1-max(abs(eig(M\N))));

% Precompute f
f = (R + Deltap) \ R(:,end);
g = f ./ (1 - f(end));

% S = Su * Sv', Sv = en
Su = R(:,end) + W(:,end) + Delta(:,end);

% System to solver A \ r
x0 = M \ full(r);
x = x0;

j = 0;

while j < maxsteps
    j = j + 1;

    xold = x;
    % x = M \ (N * x) + x0;

    % Compute l = N * x
    l = -W * x;
    % l = l - Delta * x - gamma * x;
    l = l - Delta * x + Deltap * x;
    l = l + Su * x(end);

    % Solve the linear system M*x = l
    % x = (R - gamma * eye(prod(n))) \ l; 
    x = (R + Deltap) \ l; 
    x = x - g * x(end);    

    % Update the iterate
    x = x + x0;

    if debug
        fprintf('Step %d, rel. change %e, m = %e\n', ...
            j, norm(xold - x, 1), -dot(full(pi0), x));
    end

    if norm(xold - x, 1) < tol
        break;
    end
end

m = -dot(full(pi0), x);

t = toc(timer);
y = x;
time = t;

end

