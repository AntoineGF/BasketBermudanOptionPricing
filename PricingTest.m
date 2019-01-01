%% Bermudan Basket
clear all, clc
rng default

%% Parameters
S0= [100 100 100];
K = 100;
r=0.05;
y = [0.1 0.1 0.1]; 
T=3; 
steps=30;
nsims=1000; 
sigma= [0.2 0.2 0.2];           
corr_matrix = [1 -0.25 0.25;-0.25 1	0.3;0.25 0.3 1];
nexercise = 30; 
%dummy for regression with higher power (=1 if higher powers, =0 if status quo)
high_pow=0;

%% Normal Pricing Function 
antithetic=0;
V = BasketBermudanCall(S0, K, r, y, T, steps, nsims, sigma, corr_matrix, nexercise,antithetic, high_pow);

%% Antithetic Pricing Function
antithetic=0;
V1 = BasketBermudanCall(S0, K, r, y, T, steps, nsims, sigma, corr_matrix, nexercise,antithetic, high_pow);
antithetic=1;
V_anti = BasketBermudanCall(S0, K, r, y, T, steps, nsims, sigma, corr_matrix, nexercise,antithetic, high_pow);
V = (V1+V_anti)/2;

%% Confidence intervals: number of paths
tic
high_pow=0;
antithetic=0;
nsims = [10, 50, 100, 500, 1000];
nexercises = [10,15,30];
simulations = 100; %(to compute the std and then the CI)
V = zeros(simulations, length(nsims), length(nexercises));

for e = 1:length(nexercises)
    nex = nexercises(e);
    % disp(nex)
    for j = 1:length(nsims)
        n = nsims(j);
        % disp(n)
        for i = 1:simulations
            % Removed `rng default` in CorrelatedBrownian, else no
            % confidence interval, because every simulation will yield same
            % result
            % disp(i)
            if (antithetic~=1)
            V(i,j,e) = BasketBermudanCall(S0, K, r, y, T, steps, n, sigma, corr_matrix, nex, antithetic, high_pow);
            else 
            antithetic=0;
            V1(i,j,e) = BasketBermudanCall(S0, K, r, y, T, steps, n, sigma, corr_matrix, nex, antithetic, high_pow);
            antithetic=1;
            V_anti(i,j,e) = BasketBermudanCall(S0, K, r, y, T, steps, n, sigma, corr_matrix, nex, antithetic, high_pow);           
            V = (V1+V_anti)/2;
            end
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

% Final 
subplot(3,4,[3,4,7,8,11,12])
plot(nsims, middle, '--s', 'Linewidth', 2)
xlabel('Number of simulations'), ylabel('Option Value')
legend('10 Exercise Days', '15 Exercise Days', '30 Exercise Days')
title('Comparison of the mean values')

toc
