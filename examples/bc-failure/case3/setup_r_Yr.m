% We need a different r for the _Yr cases. 
v = [3; 0; 1; 0; 0];
w = cell(1, n); 
for j = 1 : n
    if j <= l
        w{j} = tt_tensor(v(1:4)); 
    else
        w{j} = tt_tensor(v);
    end
end
r = round(ktt_kron(w{:}), 1e-8);