% We need a different r for the _Yr cases. 
v = [3; 0; 1; 0; 0]; 
w = cell(n, n); 
for j = 1 : n
    for i = 1 : n
        if j <= l
            w{i,j} = tt_tensor(ones(4,1));
        else
            w{i,j} = tt_tensor(ones(5,1));
        end
    end
    
    if j <= l
        w{j,j} = tt_tensor(v(1:4)); 
    else
        w{j,j} = tt_tensor(v);
    end
end

r = round(ktt_kron(w{1, :}), 1e-8);

for j = 2 : n
    r = round(r + ktt_kron(w{j, :}), 1e-8);
end