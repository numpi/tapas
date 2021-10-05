rng(12);

n = 160;
% A = rand(n); A = A * A';
A = triu(rand(n+1)); A = A - diag(sum(A, 2)); A = A(1:end-1,1:end-1);

bs = 2;

V = randn(n, bs);
W = randn(n, bs);

ll = W' * V;
V = V / ll;

a = {};
b = {}; c = {};

k = 60;

for j = 1 : k
    Av = A * V(:,(j-1)*bs+1 : j*bs);
    Aw = A' * W(:,(j-1)*bs+1 : j*bs);
    
    a{j} = W(:,(j-1)*bs+1 : j*bs)' * Av;
    
    Av = Av - V(:,(j-1)*bs+1 : j*bs) * a{j};
    if j > 1
        Av = Av - V(:,(j-2)*bs+1:(j-1)*bs) * b{j-1};
    end
    
    Aw = Aw - W(:,(j-1)*bs+1:j*bs) * a{j}';
    if j > 1
        Aw = Aw - W(:,(j-2)*bs+1:(j-1)*bs) *  c{j-1}';
    end
    
    nn = Aw' * Av;
    
    c{j} = nn;
    b{j} = eye(bs);
    
    Aw = Aw / b{j};
    Av = Av / c{j};
    
    V = [ V, Av ];
    W = [ W, Aw ];
    
    % T = diag(a) + diag(b(1:end-1),1) + diag(c(1:end-1),-1);    
    T = tridiag(a,b,c);
    
    S = inv(T); S(1:bs,1:bs)
end

function T = tridiag(a,b,c)
    bs = length(a{1});
    sz = length(a) * bs;
    T = zeros(sz, sz);
    
    for i = 1 : length(a)
        T((i-1)*bs+1:i*bs, (i-1)*bs+1:i*bs) = a{i};
    end
    for i = 1 : length(b)-1
        T((i-1)*bs+1:i*bs, i*bs+1:(i+1)*bs) = b{i};
    end    
    for i = 1 : length(c)-1
        T(i*bs+1:(i+1)*bs, (i-1)*bs+1:i*bs) = c{i};
    end
end

% T = diag(a) + diag(b(1:end-1),1) + diag(c(1:end-1),-1);