function [first_best, second_best] = find_most_performing(S, period)
% Input:    S: is a three dimensional array
%           period_ is the considered time period
% Output:   S1: best performing
%           S2: second-best performing
%
% Measure:  Two assets are best performing if their maximum value is greater
%           than the value of the remaining asset(s)

stock = zeros(1, size(S, 3));

% Select two best performing stocks
for i = 1:size(S, 3)
    stock(i) = max(S(period,:,i));
end

first_best = stock == max(stock);
second_best = stock == max(stock(stock < max(stock)));

end