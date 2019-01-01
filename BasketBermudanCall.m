function V = BasketBermudanCall(S0, K, r, y, T, steps, nsims, sigma, corr_matrix, nexercise,antithetic, high_pow)
% Input:    S0:             initial price vector
%           K:              strike (identical for all options)
%           r:              risk-free rate
%           y:              dividend rates vector
%           T:              Horizon
%           steps:          Number of steps
%           nsims:          Number of simulated path by stock
%           sigma:          Vector of standard deviation
%           corr_matrix:    Correlation matrix between the stocks
%           nexercise:      Number of exercise possibility
%           antithetic:     Dummy =1 if with antithetic variance reduction
%           high_pow:      	Dummy =1 if regression of higher power is used
%
% Output:   V:              option price

%% Set up
% Total periods
periods = T*steps;

% Equidistant exercise possibility
ex_intervals = ceil(periods/nexercise);

% Simulate stock price
S = CorrelatedBrownian(S0,r,y,sigma,corr_matrix,T,periods,nsims,antithetic);

% Initialize final_cash_flow
final_cash_flow = zeros(nsims, periods);

%% START LOOP BACKWARD
% First iteration is outside the loop
[first_best, ~] = find_most_performing(S, periods + 1);

% Last cash-flow
final_cash_flow(:, end) = max(S(periods+1, :, first_best) - K, 0);
% temp is only used for the regression
temp_cash_flow = final_cash_flow(:, periods);

% American periods:-1:2
% for t = periods:-1:2
for t = periods:-ex_intervals:ex_intervals
    %% 1.Find two most performing assets 
    
    [first_best, second_best] = find_most_performing(S, t + 1);
    % S1 and S2 to be used in the regression (i.e. at t-1)
    S1 = S(t,:,first_best);
    S2 =  S(t,:,second_best);
    
    % Find the ones in the money
    [S1_itm, S2_itm, position] = in_the_money(S1, S2, K);

    %% 2. Regression
    
    % Regression as proposed in the exercise
    if high_pow~=1
  
    reg_matrix = [S1_itm; S1_itm.^2; S2_itm; S2_itm.^2; S1_itm.*S2_itm]'; 
    cash_flow_regression = temp_cash_flow(position);
    dependent_var = exp(-r)*cash_flow_regression;
    
    %add vector of ones for determination of  intercept 
    x = ones(size(S1_itm,2),1);
    reg_matrix = [x, reg_matrix]; 
    
    
    % Fit regression (sometimes, more predictor than observations, hence the if)
    if size(reg_matrix,1)>5
        
        % Fit regression
        ols_fit = reg_matrix\dependent_var; 
        
        % Fitted value if hold (estimated value of holdings)
        holding_value = max(reg_matrix*ols_fit,0);
        
    elseif size(reg_matrix,1)>2
        
        % Take only two predictors due to potential error message
        reg_matrix = [S1_itm; S2_itm]'; 
        
        %add vector of ones for determination of  intercept 
        x = ones(size(S1_itm,2),1);
        reg_matrix = [x, reg_matrix]; 
        ols_fit = reg_matrix\dependent_var;
        holding_value = max(reg_matrix*ols_fit,0);
        
    else 
         
        holding_value = zeros(length(dependent_var),1); 
        
    end
    
    % If dummy high_pow=1 then we choose a regression of higher power
    elseif high_pow==1  
        
    reg_matrix = [S1_itm; S1_itm.^2; S1_itm.^3; S2_itm; S2_itm.^2; S2_itm.^3; S1_itm.*S2_itm; S1_itm.^2.*S2_itm; S1_itm.*S2_itm.^2]';  
    cash_flow_regression = temp_cash_flow(position);
    dependent_var = exp(-r)*cash_flow_regression;
    
    %add vector of ones for determination of  intercept 
    x = ones(size(S1_itm,2),1);
    reg_matrix = [x, reg_matrix]; 

    % Fit regression (sometimes, more predictor than observations, hence the if)
    if size(reg_matrix,1)>9 
        
        % Fit regression
        ols_fit = reg_matrix\dependent_var;
        
        % Fitted value if hold (estimated value of holding)s
        holding_value = max(reg_matrix*ols_fit,0);
        
    elseif size(reg_matrix,1)>5
        
        reg_matrix = [S1_itm; S1_itm.^2; S2_itm; S2_itm.^2; S1_itm.*S2_itm]';
        
        %add vector of ones for determination of  intercept 
        x = ones(size(S1_itm,2),1);
        reg_matrix = [x, reg_matrix];
        
        % Fit regression
        ols_fit = reg_matrix\dependent_var;
        
        % Fitted value if hold (estimated value of holding)s
        holding_value = max(reg_matrix*ols_fit,0);
        
    elseif size(reg_matrix,1)>2
        
        % Take only two predictors if possible
        reg_matrix = [S1_itm; S2_itm]';
        
        %add vector of ones for determination of  intercept 
        x = ones(size(S1_itm,2),1);
        reg_matrix = [x, reg_matrix];
        
        ols_fit = reg_matrix\dependent_var;
        holding_value = max(reg_matrix*ols_fit,0);
        
    else 
           
        holding_value = zeros(length(dependent_var),1); 
        
    end           
    end

    
    %% 3. Compare exercise vs hold value (strategy value)
    
    early_exercise = max(S1_itm - K, 0);
    strategy_value = compare_strategy(early_exercise, holding_value);

    % Store the result in the cash flow matrix
    final_cash_flow(position, t - 1) = strategy_value;
    
    % Get the next temporary cash-flow (for the next iteration and regression)
    cf = [final_cash_flow(:, t-1), temp_cash_flow];
    for i = 1:nsims
        if cf(i,1) == 0
            cf(i,1) = cf(i,2);
        end
    end
    temp_cash_flow = cf(:,1);
end

%% Clean the final cash flow 

% Only the first element to appear is non-zero, otherwise zero 
cleaned_cash_flow = clean_cash_flow(final_cash_flow);

% Compute the price of the option by discounting the cash-flows
discount_factor = diag(exp(-r*(1:periods)./steps));
discounted_cf = cleaned_cash_flow*discount_factor;

%% Output
V = sum(sum(discounted_cf))/nsims;
% disp(V)
