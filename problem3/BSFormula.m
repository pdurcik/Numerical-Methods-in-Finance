function V = BSFormula(t, T, S, K, r, sigma, is_call)
    % S - stock price
    % K - strike price
    % is_call - if 1, it calculates for call option, if it is not 1 it 
    %           calculates for put option
  
  
    d1 = (1./(sigma*(T-t).^(1/2))).*(log(S/K)+(r+(sigma^2/2))*(T-t));
    d2 = d1 - sigma*(T-t).^(1/2);
    PV = K*exp(-r*(T-t));
    
    % call option
    if is_call == 1
      N1c = (1/2)*erfc(-d1/(sqrt(2)));
      N2c = (1/2)*erfc(-d2/(sqrt(2)));
      V = max(N1c.*S - N2c.*PV,0); % values of call option
      V = fliplr(V);
      
    % put option (using put–call parity)
    else
      N1p = (1/2)*erfc(d1/(sqrt(2)));
      N2p = (1/2)*erfc(d2/(sqrt(2)));
      V = N2p.*PV - N1p.*S; % values of put option
      V = fliplr(V);
    end
    
end
