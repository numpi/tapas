function [ en ] = createen(m,K)
%CREATEEN Create the vector [ 0 ; ... ; 0 ; 1] in the TT format. 

if exist('K', 'var')

    en = tt_tensor([ 0 ; 1 ]);
    x = en; 

    for i = 2 : K
        en = tkron(en, x);
    end

else
    en = tt_tensor([ zeros(m(1) - 1, 1) ; 1 ]);
    
    for i = 2 : length(m)
        en = tkron(en, tt_tensor([ zeros(m(i) - 1, 1) ; 1 ]));
    end
end

end

