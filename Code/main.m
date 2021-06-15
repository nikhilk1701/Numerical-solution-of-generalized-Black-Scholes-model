clear all;
clc

r = @(S,t) 0.04;
sigma = @(S,t) 0.04;
D = @(S,t) 0.02;

T = 1;
Smax = 8;
K = 1;
h = 0.01;
k = 0.01;
S_ = 0:k:Smax;
T_ = 0:h:T;
N = T/h + 1;
M = Smax/k + 1;

U = solveEuropeanVanilla(Smax,T,K,h,k,r,D,sigma);
u_actual = zeros(N,M);

for n = 1:N
    for m = 1:M
        u_actual(n,m) = getEuropeanCallValue(S_(m),T_(n),K,r(0,0),D(0,0),sigma(0,0));
    end
end
u_actual = flipud(u_actual);
% 
figure(1)
subplot(1, 2, 1)
surf(S_,T_,U);
xlabel('S')
ylabel('t')
zlabel('V(S,t)')
title('Approximated V(S,t)')
% 
figure(1)
subplot(1,2,2)
surf(S_,T_,u_actual);
xlabel('S')
ylabel('t')
zlabel('V(S,t)')
title('Analytical Solution')
sgtitle('HODIE')

figure(2)
plot(S_, U(1,:))
xlabel('S_')
ylabel('V(S,T)')
title('Numerical approximation of value of the option at t = 0')

%% error analysis

e_max = getMaxError(U,S_,T_,K,r,D,sigma);
e_rms = getMeanSquaredError(U,S_,T_,K,r,D,sigma);
% 
figure(3)
subplot(1,2,1)
plot(T_,e_max);
xlabel('Time');
ylabel('max error at time t');
title('max error vs time t');

figure(3)
subplot(1,2,2)
plot(T_,e_rms);
xlabel('Time');
ylabel('rms error at time t');
title('rms error vs time t');
%% convergence
N = 40;
M = 64;

erms_ = zeros(4,1);
emax_ = zeros(4,1);
for i = 1:4
    h = T/N;
    k = Smax/M;
    U = solveEuropeanVanilla(Smax,T,K,h,k,r,D,sigma);
    
    S_ = 0:k:Smax;
    T_ = 0:h:T;
    erms = getMeanSquaredError(U,S_,T_,K,r,D,sigma);
    erms_(i) = erms(N);
    
    emax = getMaxError(U,S_,T_,K,r,D,sigma);
    emax_(i) = emax(N);
    
    N = 2*N;
    M = 2*M;
end

figure(4)
plot(erms_)
title('rms error vs size of mesh');

figure(5)
plot(emax_)
title('max error vs size of mesh');

prms_ = zeros(1,3);
pmax_ = zeros(1,3);

for i = 1:3
    prms_(i) = log2(erms_(i)/erms_(i+1));
    pmax_(i) = log2(emax_(i)/emax_(i+1));
end

fprintf("The convergence ratios for E_rms \n");
disp(prms_);
fprintf("The convergence ratios for E_max \n");
disp(pmax_);
