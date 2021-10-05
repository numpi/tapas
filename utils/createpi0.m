function [ pi0 ] = createpi0(m, K)
%CREATEPI0 Create the vector [ 1 ; 0 ; ... ; 0 ] in the TT format. 

if exist('K', 'var')
    pi0 = tt_tensor([ 1 ; 0 ]);
    x = pi0; 

    for i = 2 : K
        pi0 = tkron(pi0, x);
    end
else
    pi0 = tt_tensor(eye(m(1), 1));
    
    for i = 2 : length(m)
        pi0 = tkron(pi0, tt_tensor(eye(m(i), 1)));
    end
end


end

