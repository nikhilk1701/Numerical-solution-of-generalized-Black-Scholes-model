function y = getMeanSquaredError(U,S_,T,K,r_,D_,sigma_)
    r = r_(0,0);
    D = D_(0,0);
    sigma = sigma_(0,0);
    
    [N,M] = size(U);
    e_rms = zeros(N,1);
    
    for n = 1:N
        for m = 1:M
            u = getEuropeanCallValue(S_(m),T(N-n+1),K,r,D,sigma);
            e_rms(n) = e_rms(n) + (U(n,m) - u)^2;
        end
        e_rms(n) = (e_rms(n)/(M+1))^0.5;
    end
    
    y = e_rms;
end