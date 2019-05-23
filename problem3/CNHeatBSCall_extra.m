function [S,t_space, Csol] = CNHeatBSCall_extra(Sa,Sb,Numx,dt,K, r, sig, t, T)
    % calculate for biger interval (because of x->inf) nad then
    % restrict it to smaller interval (0, 3*K)
    
    % dx - spacial step
    % dt - time step
    % S - stock price discritisations
    % Sa - lowest price in scheme
    % Sb - highest price in scheme
    % K - strike price
    % Numx - number of space steps
    % X - spacial discritisation
    % t_space - time discritisation
    % Usol - matrix of the values of u(x,t)
    
    if Sa >0            % If Sa is positive 
        a = log(Sa/K);  % the lowest spacial value is defined as log(Sa) because x=ln(S/K)
    else                % otherwise we define
        a = -5;         % the lowest spacial value is defined as -5
    end
    be = log(Sb/K)*1.5; % the highest spacial value is defined as log(Sb) because x=ln(S/K)
    b = log(Sb/K);      % the highest spacial value is defined as log(Sb) because x=ln(S/K)
    dx = (b-a)/Numx;    % the spacial step is defined as the distance between 
                        % a and b divided by the number of steps Numx
    
    tauI = 0;               % tauI is the lowest value of tau and is defined as 0
    tauF = (sig^2)*(T-t)/2; % tauF is the highest value of tau and is defined as (sig^2)*(T-t)/2
    dtau = dt*(sig^2)/2;    % the temporal step in the variable tau
    
    X = (a:dx:be)';         % define X as a spacial vector between a and b with spacial steps dx
                            % i.e: X = |a  a+dx  a+2*dx  a+3*dx  ......  a+(M-2)*dx  b|'
    S = K*exp((a:dx:b)');   % S is defined as the price vector with stock prices between Sa and Sb
                            % given by S = |Sa  K*exp(a+dx)  K*exp(a+2*dx)  ......  K*exp(a+(M-2)*dx)  Sb|'
    tau = tauI:dtau:tauF;   % tau is defined as a time vector between tauI and tauF with time steps dtau
                            % i.e: tau = |0  dt  2*dt  3*dt  ...... (N-2)*dt  T|
    t_space = t:dt:T;       % t_space is defined as a time vector between t and T with time steps dt
                            % i.e: t_space = |t  t+dt  t+2*dt  t+3*dt  ...... t+(N-2)*dt  T|
                            
    lambda = dtau/(dx^2);       % define lambda as dtau/(dx^2) 
    M = length(X);              % define M as the length of vector X
    N = length(tau);            % define N as the length of vector tau
    
    k1 = 2*r/(sig^2);           % k1 is defined as 2*r/(sig^2) 
    
   
    A = tridiag(M,1-lambda, lambda/2, lambda/2);    % A is defined as a M dimensional tridiagonal matrix
                                                    % with central diagonal values of 1-lambda
                                                    % upper and lower diagonal values of lambda/2
                                                    
    % A is modified in order for the boundary conditions of u to not be affected
    % the first row of A is defined as |1 0 0 0 ..... 0 0|
    A(1,1) = 1; A(1,2) = 0;
    % the last row of A is defined as |0 0 0 0 ..... 0 1|
    A(end, end-1) = 0; A(end, end) = 1;
                                                    
    % Hence A is a MxM matrix of the form
    %
    %     |     1           0          0           0       .........      0           0           0     |
    %     | lambda/2    1-lambda    lambda/2       0       .........      0           0           0     |
    %     |     0       lambda/2    1-lambda    lambda/2   .........      0           0           0     |
    %     |     0           0       lambda/2    1-lambda   .........      0           0           0     |
    % A = |     .           .          .           .       .........      .           .           .     |
    %     |     .           .          .           .       .........      .           .           .     |
    %     |     .           .          .           .       .........      .           .           .     |
    %     |     0           0          0           0       .........   1-lambda    lambda/2       0     |
    %     |     0           0          0           0       .........   lambda/2    1-lambda    lambda/2 |
    %     |     0           0          0           0       .........      0           0           1     |
    

    
    B = tridiag(M,1+lambda, -lambda/2, -lambda/2);  % B is defined as a M dimensional tridiagonal matrix
                                                    % with central diagonal values of 1+lambda
                                                    % upper and lower diagonal values of -lambda/2   
    
    % B is modified in order for the boundary conditions of u to not be affected
    % the first row of B is defined as |1 0 0 0 ..... 0 0|
    B(1,1) = 1; B(1,2) = 0;
    % the last row of B is defined as |0 0 0 0 ..... 0 1|
    B(end, end-1) = 0; B(end, end) = 1;
  
    % Hence B is a MxM matrix of the form
    %
    %     |     1           0          0           0       .........      0           0           0     |
    %     |-lambda/2    1+lambda   -lambda/2       0       .........      0           0           0     |
    %     |     0      -lambda/2    1+lambda   -lambda/2   .........      0           0           0     |
    %     |     0           0      -lambda/2    1+lambda   .........      0           0           0     |
    % B = |     .           .          .           .       .........      .           .           .     |
    %     |     .           .          .           .       .........      .           .           .     |
    %     |     .           .          .           .       .........      .           .           .     |
    %     |     0           0          0           0       .........   1+lambda   -lambda/2       0     |
    %     |     0           0          0           0       .........  -lambda/2    1+lambda   -lambda/2 |
    %     |     0           0          0           0       .........      0           0           1     |
    
        
    Usol = zeros(M,N);  % the matrix of solutions is originaly defined as a
                        % NxM matrix of only zeroes
    
    Usol(:,1) = max(exp(((k1+1)*X)/2) - exp(((k1-1)*X)/2), 0);  % the values of u(x,tau) when tau = 0 are defined as
                                                                % u(x,0) = max(exp(((k1+1)*X)/2) - exp(((k1-1)*X)/2), 0)
    
    for i=2:N
        Usol(:,i) = B\(A*Usol(:,i-1));  % Usol(:,n) is defined as (B^(-1))*A*Usol(:,n-1)
    end                                 % because  B*Usol(:,n) = A*U(:,n-1) for any n
    
    
    [x_val, tau_val] = meshgrid(X,tau);
    
    C1 = K*(exp((-(k1-1)*x_val')/2 -((k1+1)^2)*tau_val'/4)).*Usol;    % To represent the change of variables C is
                                                            % defined as K*(exp((-(k1-1)*X)/2 -((k1+1)^2)*tau/4)).*Usol
    
    C2 = fliplr(C1);  % the matrix C is horizontaly fliped because 
                      % t and tau directions are opposed to one another
    
    n = length(S);
    Csol = C2(1:n,:);
                    

    % Therefore Csol a matrix is of the form:
    %
    %        |      C(Sa,t)                 C(Sa,t+1*dt)                 C(Sa,t+2*dt)                 C(Sa,t+3*dt)         .........         C(Sa,T)       |
    %        | C(K*exp(a+1*dx),t)      C(K*exp(a+1*dx),t+1*dt)      C(K*exp(a+1*dx),t+2*dt)      C(K*exp(a+1*dx),t+3*dt)   .........    C(K*exp(a+1*dx),T) |
    %        | C(K*exp(a+2*dx),t)      C(K*exp(a+2*dx),t+1*dt)      C(K*exp(a+2*dx),t+2*dt)      C(K*exp(a+2*dx),t+3*dt)   .........    C(K*exp(a+2*dx),T) |
    %        | C(K*exp(a+3*dx),t)      C(K*exp(a+3*dx),t+1*dt)      C(K*exp(a+3*dx),t+2*dt)      C(K*exp(a+3*dx),t+3*dt)   .........    C(K*exp(a+3*dx),T) |
    % Csol = |         .                         .                            .                            .               .........            .          |
    %        |         .                         .                            .                            .               .........            .          |
    %        |         .                         .                            .                            .               .........            .          |
    %        |         .                         .                            .                            .               .........            .          |
    %        |      C(Sb,t)                C(Sb,t+1*dt)                  C(Sb,t+2*dt)                 C(Sb,t+3*dt)         .........         C(Sb,T)       |
                    
end




% tridiag(n, c, l, u) creates a triagonal matrix of size n with central
% values c lower lower values l and upper values u

function M = tridiag(n, c, l, u)
    nOnes = ones(n, 1);                                                         % create a n length vector of 1's
    
    M = diag(c*nOnes) + diag(l*nOnes(1:n-1), -1) + diag(u*nOnes(1:n-1), 1);     % create a diagonal matrix with value c in the diagonal
end                                                                             % a upper diagonal matrix with value d and a lower diagonal
                                                                                % matrix with value l and sum the three matrices
