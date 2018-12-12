function cleaned_cash_flow = clean_cash_flow(final_cash_flow)
% Input     final_cash_flow: the cash flow matrix to be cleand (only one 
%           non zero element by row
% output    final_cash_flow: is cleaned

% Iterate to eliminate the unwelcomed elements
for i = 1:size(final_cash_flow,1)
    temp = final_cash_flow(i, :);
    pos = find(temp ~= 0, 1);
    final_cash_flow(i, pos+1:end) = 0;
end

cleaned_cash_flow = final_cash_flow;