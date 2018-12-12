function [S1_itm, S2_itm, position] = in_the_money(S1, S2, K)
% Input     S1, S2: two best performing assets 
%           K: strike price
% Output    S1_itm, S2_itm: the two stock prices in the money
%           position: the index position of the price in-the-money
% Criterion We take the same index for S2 as the ones that are itm in S1

% Select prices that are in the money (for best performing one)
S1_itm = S1(max(S1 - K, 0) > 0);
S2_itm = S2(max(S1 - K, 0) > 0);
position = find(max(S1 - K, 0) > 0);
end

