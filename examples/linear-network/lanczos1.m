rng(12);

n = 1600;
% A = rand(n); A = A * A';
A = triu(rand(n+1)); A = A - diag(sum(A, 2)); A = A(1:end-1,1:end-1);

V = randn(n,1);
W = randn(n,1);

ll = dot(V, W);
V = V; W = W / ll;

a = [];
b = []; c = [];

k = 120;

for j = 1 : k
    Av = A * V(:,j);
    Aw = A' * W(:,j);
    
    a(j) = dot(W(:,j), Av);
    
    Av = Av - a(j) * V(:,j);
    if j > 1
        Av = Av - b(j-1) * V(:,j-1);
    end
    
    Aw = Aw - a(j) * W(:,j);
    if j > 1
        Aw = Aw - c(j-1) * W(:,j-1);
    end
    
    nn = dot(Aw, Av);
    
    c(j) = sqrt(abs(nn));
    b(j) = nn / c(j);
    
    Aw = Aw / b(j);
    Av = Av / c(j);
    
    V = [ V, Av ];
    W = [ W, Aw ];
    
    T = diag(a) + diag(b(1:end-1),1) + diag(c(1:end-1),-1);    
    
    S = inv(T); S(1,1)
end

T = diag(a) + diag(b(1:end-1),1) + diag(c(1:end-1),-1);