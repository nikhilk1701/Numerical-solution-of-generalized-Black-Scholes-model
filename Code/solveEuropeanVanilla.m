function y = solveEuropeanVanilla(Smax,T,K,h,k,r,D,sigma)
    S_ = 0:k:Smax;
    M = Smax/k + 1;
    N = T/h + 1;
    T = 0:h:T;
    U = zeros(N,M);
    
    a0 = @(n,m) -r(S_(m),T(n));
    a1 = @(n,m) (r(S_(m),T(n)) - D(S_(m),T(n)))*S_(m);
    a2 = @(n,m) (sigma(S_(m),T(n))^2*S_(m)^2)/2;
    f = @(n,m) 0;
    
    % todo properly fill g1 and g2
    g1 = @(t) 0;
    g2 = @(t) Smax*exp(-D(Smax,t)*t) - K*exp(-r(Smax,t)*t);
    
    eps = 1e-6;
    for i = 1:M
        %U(1,i) = max(0, (i-1)*k - K);
        if (i-1)*k - K <= -eps
            U(1, i) = 0;
        elseif -eps < (i-1)*k - K && (i-1)*k - K < eps
            temp = (i-1)*k - K;
            temp_1 = temp/eps;
            U(1, i) = (35*eps/256) + (1/2)*temp + (35/64)*temp*temp_1 + (-35/128)*temp*temp_1^3 + (7/64)*temp*temp_1^5 + (-5/256)*temp*temp_1^7;
        else
            U(1, i) = (i-1)*k - K;
        end
    end
    %disp(U(1,2));
    
    beta1 = @(n,m) (6*h*a2(n,m+1) + 2*h^2*a1(n,m+1)) / (6*h*a2(n,m+1) + 2*h^2*a1(n,m+1) + h^2*a1(n,m));
    beta2 = @(n,m) 1 - beta1(n,m);
    alpha_m = @(n,m) (beta1(n,m)*(-2*a2(n,m) + h*a1(n,m)) + beta2(n,m)*(-2*a2(n,m+1) - h*a1(n,m+1))) / (2*h^2);
    alpha_p = @(n,m) (beta1(n,m)*(-2*a2(n,m) - h*a1(n,m)) + beta2(n,m)*(-2*a2(n,m+1) - 3*h*a1(n,m+1) - 2*h^2*a0(n,m+1))) / (2*h^2);
    alpha_c = @(n,m) (beta1(n,m)*(4*a2(n,m) - 2*h^2*a0(n,m)) + beta2(n,m)*(4*a2(n,m+1) + 4*h*a1(n,m+1))) / (2*h^2);
    
    % todo for n = 2
    A = zeros(M,M);
    b = zeros(M,1);
    n = 2;
    for m = 1:M
        if m == 1
            A(m,m) = 1;
            b(m) = g1(T(n));
        elseif m < M
            A(m,m-1) = alpha_m(n,m);
            A(m,m) = alpha_c(n,m) + beta1(n,m)/k;
            A(m,m+1) = alpha_p(n,m) + beta2(n,m)/k;
            
            b(m) = beta1(n,m)*(f(n,m) + U(n-1,m)/k);
            b(m) = b(m) + beta2(n,m)*(f(n,m+1) + U(n-1,m+1)/k);
        else
           A(m,m) = 1;
           b(m) = g2(T(n));
        end
    end
    %disp(g2(T(n)));
    U(n,:) = (A\b)';
    %disp(U(n,:));
    
    for n = 3:N
        A = zeros(M,M);
        b = zeros(M,1);
        for m = 1:M
            if m == 1
                A(m,m) = 1;
                b(m) = g1(T(n));
            elseif m < M
                A(m,m-1) = alpha_m(n,m);
                A(m,m) = alpha_c(n,m) + 1.5*beta1(n,m)/k;
                A(m,m+1) = alpha_p(n,m) + 1.5*beta2(n,m)/k; 
                
                b(m) = beta1(n,m)*(f(n,m) + 2*U(n-1,m)/k - 0.5*U(n-2,m)/k);
                b(m) = b(m) + beta2(n,m)*(f(n,m+1) + 2*U(n-1,m+1)/k - 0.5*U(n-2,m+1)/k);
            else
                A(m,m) = 1;
                b(m) = g2(T(n));
            end
        end
        U(n,:) = (A\b)';
        %disp(g2(T(n)));
    end
    y = flipud(U);
end