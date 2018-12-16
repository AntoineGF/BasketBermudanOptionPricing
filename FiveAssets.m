%% Bermudean Basket
clear all, clc

%% 10 assets Parameters
S0= ones(1,5)*100;
K = 100;
r=0.05;
y = ones(1,5)*0.1;
T=3; 
steps=30;
nsims=100; 
sigma= [0.2 0.1 0.2 0.1 0.3];
nexercise = 30; 
corr_matrix = [
    1 0.2 -0.1 0.1 0;
    0.2 1 0.3 0.2 -0.1;
    -0.1 0.3 1 -0.3 0.25; 
    0.1 0.2 -0.1 1 0;
    0 -0.1 0.25 0 1;];
%% Function 
V = BasketBermudanCall(S0, K, r, y, T, steps, nsims, sigma, corr_matrix, nexercise);
disp(V)
%% Confidence intervals: number of paths
tic
nsims = [10, 50, 100];
nexercises = [10,15,30];
simulations = 50; %(to compute the std and then the CI)
V = zeros(simulations, length(nsims), length(nexercises));

for e = 1:length(nexercises)
    nex = nexercises(e);
    % disp(nex)
    for j = 1:length(nsims)
        n = nsims(j);
        % disp(n)
        for i = 1:simulations
            % Removed `rng default` in CorrelatedBrownian
            % disp(i)
            V(i,j,e) = BasketBermudanCall(S0, K, r, y, T, steps, n, sigma, corr_matrix, nex);
        end
    end
end

% Confidence interval 
upper = zeros(length(nsims), length(nexercises));
lower = upper;
middle = upper;
for i = 1:length(nexercises)
    middle(:,i) = mean(V(:,:,i));
end

% Compute the bounds
for j = 1:length(nexercises)
    for i = 1:length(nsims)
        upper(i, j) = middle(i,j) + 1.96*std(V(:,i,j))/sqrt(nsims(i));
        lower(i, j) = middle(i,j) - 1.96*std(V(:,i,j))/sqrt(nsims(i));
    end
end

% Plots (only for 3 cases) (fancy idx, ix are just for the subplot)
plot_titles = ["10 Exercise Days", "15 Exercise Days", "30 Exercise Days"];
idx = [0, 3, 6];
for i = 1:3
    ix = [i + idx(i), i + 1 + idx(i)];
    subplot(3,4, ix)
    plot(nsims, middle(:,i), '--rs', 'Linewidth', 2), grid on, hold on
    plot(nsims, upper(:,i), 'b-')
    plot(nsims, lower(:,i), 'b-')
    legend('Mean', 'Upper bound', 'Lower bound')
    xlabel('Number of simulations'), ylabel('Option Value')
    title(plot_titles(i))
end

% Finally 
subplot(3,4,[3,4,7,8,11,12])
plot(nsims, middle, '--s', 'Linewidth', 2)
xlabel('Number of simulations'), ylabel('Option Value')
legend('10 Exercise Days', '15 Exercise Days', '30 Exercise Days')
title('Comparison of the mean values')

% The less the power the option give, the cheaper it is (tends to european
% price, whereas more exercie possibilities make the price tend to an
% american option)
toc
