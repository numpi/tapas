function [m, time, y] = inv_tt_expsumst(pi0, r, R, W, absorbing_states, shift, ttol, tol)
%INV_TT_EXPSUMST 

% MTTD Q expsums tt
timer = tic;

debug = true;
maxsteps = inf;
interval_report = 1;

[~, Delta] = ktt_infgen(R, W, ttol);
[S, DeltapC, Deltap, Wsync] = inv_computeS(R, W, Delta, absorbing_states, shift, ttol);

A2 = round(Wsync + (Delta - Deltap), ttol);

k = length(R);

DAt = R;
for j = 1 : k
    DAt{j} = (- DAt{j} + DeltapC{j})';
end

% Scale everything
scl = minmaxeig(DAt);

for j = 1 : k
    DAt{j} = DAt{j} / scl;
end
S = S / scl;
A2 = A2 / scl;

expn = ceil(- 6 * log(ttol)/pi);


% r = ktt_ones(n) - en;
y0 = pi0;
m = ttexpsummldivide(DAt, y0, expn, ttol);
y = m;
nrmY = norm(y);
j = 1;
rho = 1;
while j < maxsteps
    y = ttexpsummldivide(DAt, (A2 - S)' * y, expn, ttol);
    m = round(m + y, ttol);
    oldnrmY = nrmY;
    nrmY = norm(y); nrmM = norm(m);
    oldrho = rho;
    rho = nrmY / oldnrmY;    
    err = nrmY / (1 - rho) / nrmM;
    if debug && mod(j, interval_report) == 0
        fprintf('Step %d, Neumann residue ~ %e, norm(m) = %e, rank = %f, rank y = %f, spectral radius ~ %e\n', ...
            j, nrmY, nrmM, max(rank(m)), max(rank(y)), rho);
        fprintf('Measure estimate: %e (err. estimate = %e, est. upper bound = %e)\n', ...
            dot(m, r) / scl, err, dot(m, r) * (1 + err) / scl);
    end

    if rho < 1 && rho <= oldrho * (1 + 1e-2) && err < tol
        break
    end

    j = j + 1;
end
y = y / scl;
m = m / scl;
t = toc(timer);
time = t;

m = dot(r, m);

end

