function y = getMaxError(U,S_,T,K,r_,D_,sigma_)
    r = r_(0,0);
    D = D_(0,0);
    sigma = sigma_(0,0);
    
    [N,M] = size(U);
    e_max = zeros(N,1);
    
    for n = 1:N
        for m = 1:M
            u = getEuropeanCallValue(S_(m),T(N-n+1),K,r,D,sigma);
            e_max(n) = max(e_max(n), abs(U(n,m) - u));
        end
    end
    
    y = e_max;
end