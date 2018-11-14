%% Function to compute the price of an american put 
% Longstaff/Schwartz method

%% Generate Price paths 
clear all, clc
% For Reproducibility
rng default;

% Parameters (not same as exercise)
T = 1;          % Time to maturiy
N = 3;         % Number of periods
path = 8;        % Number of paths 
sigma = 0.1;    % Volatility
r = 0.06;       % Risk-free return
y = 0;          % Dividend yield
S0 = 1;       % Price at t0 (same for all assets)
M = 1;          % Number of stocks
K = 1.1;

% Check the option price response to sigma and r increases
sigma = linspace(0.1, 0.4, 25);
r = linspace(0, 0.1, 25); 
K = linspace(1.1, 1.4, 25);
price = zeros(length(sigma),3);
for i = 1:length(sigma)
    price(i,1) = SingleAmericanPut(T, N, path, sigma(i), 0.06, y, S0, 1.1);
    price(i,2) = SingleAmericanPut(T, N, path, 0.2, r(i), y, S0, 1.1);
    price(i,3) = SingleAmericanPut(T, N, path, 0.2, 0.06, y, S0, K(i)); 
end

plot(price)
legend('Increase Sigma', 'Increase Drift', 'Increase K')

% Different strike, different volatilies
sigma = linspace(0.1, 0.4, 25);
price = zeros(length(sigma),3);
for i = 1:length(sigma)
    price(i,1) = SingleAmericanPut(T, N, path, sigma(i), 0.06, y, S0, 1.1);
    price(i,2) = SingleAmericanPut(T, N, path, sigma(i), 0.06, y, S0, 1.2);
    price(i,3) = SingleAmericanPut(T, N, path, sigma(i), 0.06, y, S0, 1.4); 
end

plot(sigma, price), grid on
xlim([0.1,0.4])
xlabel('Volatility')
ylabel('Option price')
legend('K = 1.1', 'K = 1.2', 'K = 1.4')
title('American Put Option')





%% PREPARATION CODE
x = zeros(n+1, path);
for j = 1:path
    x(:,j) = geometric_brownian(n, r, sigma, T);
end

stock_paths = x';
plot(x)
 
% Initialize some of the matrices (value -> strategy values; type -> 1 for
% hold, 0 else).
strategy_value = zeros(path,n);
strategy_type = zeros(path,n); % 1 for hold, 0 for else; last column always empty
% To stack the results to then do the regression
stacked = zeros(path, 1);

stock_paths = [1	1.09	1.08	1.34
1	1.16	1.26	1.54
1	1.22	1.07	1.03
1	0.93	0.97	0.92
1	1.11	1.56	1.52
1	0.76	0.77	0.9
1	0.92	0.84	1.01
1	0.88	1.22	1.34];

%% While loop to go backwards
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

%% Calculate the payoffs
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
disp(V0)
