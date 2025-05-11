clc; clear; close all;

%% 1) Setup
n = 10; h = 24;
max_shift = 0.5;
gamma     = 0.01;
lambda    = 0.3;  % <- optimized
alpha_min = 1;
base_price= 5;

% Baseline load shape
base_curve = [0.3 0.3 0.3 0.3 0.3 0.4 0.6 0.8 0.7 0.6 0.5 0.5 ...
              0.5 0.5 0.6 0.8 1.0 1.2 1.5 1.6 1.4 1.2 0.8 0.5];
rng(0);
uf = 0.9 + 0.2*rand(n,1);
E0 = repmat(base_curve,n,1).*repmat(uf,1,h);
L_base = sum(E0,1);

% Prices and incentives
rho = 6*ones(1,h);
rho(1:6)   = 5;
rho(19:22) = 9;
alpha = zeros(n,h);
alpha(:,19:22) = 1.5;
pct_dev = (rho - base_price) ./ base_price;

%% 2) Elasticity matrices (optimized)
Ep = -0.3*eye(h) + 0.1*(ones(h)-eye(h));  % <- optimized
Ea =  0.8*eye(h) - 0.05*(ones(h)-eye(h)); % <- optimized
Ecomb = lambda*Ep + (1-lambda)*Ea;        % combined elasticity

resp_vec = Ecomb * pct_dev.';             % h x 1

%% 3) Hybrid Model: Load Shift Calculation with Relaxed Neutrality
DE3 = zeros(n,h);
net_shift_limit_frac = 0.02;  % allow 2% net shift per user
for i = 1:n
    for j = 1:h
        DE3(i,j) = E0(i,j) * resp_vec(j);
    end
    % Cap per hour
    DE3(i,:) = max(min(DE3(i,:), max_shift * E0(i,:)), -max_shift * E0(i,:));
    % Relaxed energy neutrality (±2% of total load)
    total_change = sum(DE3(i,:));
    net_limit = net_shift_limit_frac * sum(E0(i,:));
    if abs(total_change) > net_limit
        DE3(i,:) = DE3(i,:) - total_change / h; % redistribute
    end
end

L3 = sum(E0 - DE3,1);
profit3 = sum(sum((repmat(rho,n,1) - alpha).*DE3)) - gamma*sum(sum(DE3.^2));
DR3     = sum(DE3,1);

%% 4) Plot Aggregate Load Profiles with Shift Highlight
figure;
plot(1:h, L_base, 'k-', 'LineWidth', 2); hold on;
plot(1:h, L3, 'r--', 'LineWidth', 2);

h1 = area(1:h, max(0, L_base - L3)); 
set(h1, 'FaceColor', [0.9 0.7 0.7], 'EdgeColor', 'none');

h2 = area(1:h, max(0, L3 - L_base)); 
set(h2, 'FaceColor', [0.7 0.9 0.7], 'EdgeColor', 'none');

legend('Baseline','Hybrid','Reduced','Increased','Location','best');
xlabel('Hour','FontWeight', 'bold', 'FontSize', 12); ylabel('kWh','FontWeight', 'bold', 'FontSize', 12);
title('Hybrid DR: Load Reduction and Increase','FontWeight', 'bold', 'FontSize', 12);
grid on;

%% 5) Profit and Demand Shift Stats
peak_hours = 19:22;
offpeak_hours = [1:6 23 24];
DR_peak = sum(sum(DE3(:,peak_hours)));
DR_off  = sum(sum(DE3(:,offpeak_hours)));

fprintf('\n--- HYBRID OPTIMIZED RESULTS ---\n');
fprintf('Total ESP Profit       : Rs %.2f\n', profit3);
fprintf('Peak Hour Reduction    : %.2f kW\n', abs(DR_peak));
fprintf('Off-Peak Hour Increase : %.2f kW\n', DR_off);

%%%% 6) Optional: Plot Net Shifting

% 6a) Compute per-user net energy shift
user_net_shift = sum(DE3, 2); % sum over hours per user

% 6b) Compute user-level profits
user_profit = zeros(n,1);
for i = 1:n
    user_profit(i) = sum((rho - alpha(i,:)) .* DE3(i,:)) - gamma * sum(DE3(i,:).^2);
end

% 6c) Display table manually

fprintf('User\tProfit_Rs\n');
for i = 1:n
    fprintf('%2d\t%10.3f\t%10.2f\n', i, user_profit(i));
end

% 6d) Plot per-user net shift
figure;
bar(1:n, user_net_shift);
title('Net hourly Demand Shift per User','FontWeight', 'bold', 'FontSize', 12);
xlabel('User'); ylabel('Net Demand Shift per user (kW)');
grid on;

% 6e) Plot aggregate shift per hour
figure;
bar(1:h, DR3);
title('Net Aggregate Hourly Demand Shift (Hybrid)','FontWeight', 'bold', 'FontSize', 12);
xlabel('Hour','FontWeight', 'bold', 'FontSize', 12); ylabel('Net Demand Shift in(kW)','FontWeight', 'bold', 'FontSize', 12);
grid on;
%% 6f) Plot Hourly and Cumulative ESP Profit

% Hourly profit (vector of 24 values)
hourly_profit = sum((repmat(rho, n, 1) - alpha) .* DE3, 1) - gamma * sum(DE3.^2, 1);

% Cumulative profit
cumulative_profit = cumsum(hourly_profit);

% Plot
figure;

bar(1:h, hourly_profit);
title('Energy Retailer Hourly Profit','FontWeight', 'bold', 'FontSize', 12);
xlabel('Hour','FontWeight', 'bold', 'FontSize', 12); ylabel('Profit (Rs)','FontWeight', 'bold', 'FontSize', 12);
grid on;

