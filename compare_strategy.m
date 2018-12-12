function strategy_value = compare_strategy(early_exercise, holding_value)
% Input     early_exercise: early exercise values 
%           holding vlaue: fitted value from the regression
% Output    

strategy_value = zeros(length(early_exercise), 1);

for i = 1:length(early_exercise)
    if early_exercise(i) > holding_value(i)
        strategy_value(i) = early_exercise(i);
    else 
        strategy_value(i) = 0;
    end
end

end
