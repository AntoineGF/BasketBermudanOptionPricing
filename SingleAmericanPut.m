%% Function to compute the price of an american put 
% Longstaff/Schwartz method

function V0 = SingleAmericanPut(T, N, path, sigma, r, y, S0, K)
% INPUTS:
% T:        Time to maturiy
% N:        Number of periods
% path:     Number of paths 
% sigma:    Volatility
% r:        Risk-free return (drift)
% y:        Dividend yield
% S0:       Price at t0 (same for all assets)
% K:        Strike price

% OUTPUT: 
% V0:       Option price

% Generate Stock Path
rng default;
x = zeros(N+1, path);
for j = 1:path
    x(:,j) = geometric_brownian(N, r, sigma, T);
end

% plot
stock_paths = x';
plot(x), grid on
title('Stock Paths')

% Initialize some of the matrices (value -> strategy values; type -> 1 for
% hold, 0 else).
strategy_value = zeros(path,N);
strategy_type = zeros(path,N); % 1 for hold, 0 for else; last column always empty
% To stack the results to then do the regression
stacked = zeros(path, 1);

% While loop to go backwards
% Count the iterations 
n = N;
i = 1; 
while n > 1
    
    % Have to distinguish between first iteration and subsequent ones
    if i == 1
        % Compute payoffs of the last period
        payoffs = max(K - stock_paths(:,n+1),0);

        % Keep only in-the-money options in previous period and discount it
        % (used as Y in the regression)
        pos = find(max(K- stock_paths(:,n), 0));
        
        % Discounted payoffs
        payoffs_itm = exp(-r) * payoffs(pos); 
        
        % Exercising strategy in t-1
        payoffs_exercise = max(K- stock_paths(:,n), 0);
        
        % Regression (using stocks in previous period)
        S = stock_paths(pos,n);
        S_squared = S.^2;
        reg_matrix = [S, S_squared];
        % Fit regression
        ols_fit = regstats(payoffs_itm, reg_matrix, 'linear');
        % Fitted value if hold (estimated value of holding)s
        hold_value = ols_fit.yhat;

        % Compare the resulting values (hold vs exercise early) 

        % Here for when the holding strategy is best
        % Note: 2 conditions because it must be better than exercising 
        % AND in-the-money in the last period. 
        pos_hold = pos(hold_value > payoffs_exercise(pos));
        strategy_value(pos_hold, n) = payoffs(pos_hold);
        strategy_type(pos_hold, n-1) = ones(length(pos_hold), 1);

        % Here for when the early-exercising is best
        pos_exercise = pos(payoffs_exercise(pos) > hold_value);
        strategy_value(pos_exercise, n-1) = payoffs_exercise(pos_exercise);

        % Then, stack the two columns (in t and t-1) to continue for the next
        % iteration. 
        stacked(strategy_value(:,n) ~= 0) = strategy_value(strategy_value(:,n) ~= 0, n);
        stacked(strategy_value(:,n-1) ~= 0) = strategy_value(strategy_value(:,n-1) ~= 0, n-1);

        % stacked has to be used for the next iteration in the regression

    else 
        % Compute payoffs of the current period (exercising strategy)
        payoffs_exercise = max(K - stock_paths(:,n),0);

        % Keep only in-the-money options in previous period and discount it
        % (used as Y in the regression)
        pos = find(max(K- stock_paths(:,n), 0));
        
        % Discounted payoffs
        payoffs_itm = exp(-r) * stacked(pos);
        
        % Regression (using stocks in previous period)
        S = stock_paths(pos,n);
        S_squared = S.^2;
        reg_matrix = [S, S_squared];
        % Fit regression
        ols_fit = regstats(payoffs_itm, reg_matrix, 'linear');
        % Fitted value if hold (estimated value of holding)s
        hold_value = ols_fit.yhat;

        % Compare the resulting values (hold vs exercise early) 
        % Before: have the payoffs of the subsequent period. 
        % If not in the money, then we will NOT keep the result in the
        % stacked list.
        % This will tell you if we have to remove a holding strategy to do
        % the next regression
        % test = stacked;
        % keep_hold = stacked == 0;
        
        % Here for when the holding strategy is best
        % Note: 2 conditions because it must be better than exercising 
        % AND in-the-money in the following period. (see example)
        pos_hold = pos(hold_value > payoffs_exercise(pos));
        strategy_value(pos_hold, n-1) = stacked(pos_hold);
        strategy_type(pos_hold, n-1) = ones(length(pos_hold), 1);

        % Here for when the early-exercising is best
        pos_exercise = pos(payoffs_exercise(pos) > hold_value);
        strategy_value(pos_exercise, n-1) = payoffs_exercise(pos_exercise);

        % Then, stack the two columns (in t and t-1) to continue for the next
        % iteration. 
        stacked(strategy_value(:,n) ~= 0) = strategy_value(strategy_value(:,n) ~= 0, n);
        stacked(strategy_value(:,n-1) ~= 0) = strategy_value(strategy_value(:,n-1) ~= 0, n-1);

        % stacked has to be used for the next iteration in the regression  
        % BUT also have to check that a hold strategy will be in the-money
        % that's where we use keep_hold boolean
        % stacked = stacked(~(keep_hold & (strategy_type(:,n-1) == 1)));
        
    end
    
    % CONTROLS 
    n = n - 1; %backwards
    i = i + 1; %count iterations  
end

% Calculate the payoffs
% Scenario in which we don't keep the payoffs
% If it's an "hold" strategy AND payoff = 0 at the end
% Else, we keep. 
cash_flows = strategy_value; 
% Second loop to keep only the previous element
for i = 1:path
    for j = 2:N
        if strategy_value(i,j-1) ~= 0
            cash_flows(i, j) = 0;
        end
    end
end

% Compute the returns
V0 = zeros(n,1);
for j = 1:N
    V0(j) = sum(cash_flows(:,j)) * exp(-j*r);
end

V0 = sum(V0) / path;

end

