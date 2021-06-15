function y = getEuropeanCallValue(S,t,K,r,D,sigma)
    d1 = @(S,t) (log(S/K) + (r - D + sigma^2/2)*t) / (sigma*(t^0.5));
    d2 = @(S,t) d1(S,t) - sigma*(t^0.5);
    u = @(S,t) S*exp(-D*t)*normcdf(d1(S,t)) - K*exp(-r*t)*normcdf(d2(S,t));
    
    y = u(S,t);
end
    