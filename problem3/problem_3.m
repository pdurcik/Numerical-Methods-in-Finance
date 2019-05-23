% Run this file for the solutions of problem 3.

% Parameters:
r = 0.1;
sigma = 0.2;
K = 30;
dt = 0.2;
dx = 1;
a = 10^(-5);
b = 3*K;
t = 0;
T = 1;

% Crank-Nicolson scheme:
x_numb = 100;
[s1,t1,C1] = CNHeatBSCall_extra(a,b,x_numb,dt,K,r,sigma,t,T);

[T_grid,S_grid] = meshgrid(t1,s1);

figure;
subplot(1,2,1);
surf(T_grid,S_grid,C1)
title('Crank-Nicolson')
xlabel('Time')
ylabel('Stock price')
zlabel('Call price')

% Black-Scholes formula:
is_call = 1;
C2 = BSFormula(t,T_grid,S_grid,K,r,sigma,is_call);

subplot(1,2,2);
surf(T_grid,S_grid,C2)
title('Black-Scholes')
xlabel('Time')
ylabel('Stock price')
zlabel('Call price')

% Error function and L_inf norm:
Cerror = abs(C2-C1);

figure;
surf(T_grid,S_grid,Cerror)
title('Error function')
xlabel('Time')
ylabel('Stock price')
zlabel('Error value')

error_norm = norm(Cerror, inf);
disp(sprintf('The L_inf norm of the error function is %d.', error_norm))