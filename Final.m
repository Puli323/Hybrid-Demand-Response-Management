clc; clear; close all;

%% 1) Setup
n = 10; h = 24;
max_shift = 0.5;
gamma     = 0.01;
lambda    = 0.5;
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

%% 2) Incentive-Only Model (no elasticity matrix)
% elasticity by hour (from your image)
xi = zeros(1,h);
xi([1:6 23:24]) = 0.5;  % 00–06 & 22–23
xi(7:16)        = 0.3;  % 07–16
xi(17:22)       = 0.1;  % 17–21

DE1 = zeros(n,h);
for i = 1:n
    for j = 1:h
        DE1(i,j) = E0(i,j) * xi(j) * ((alpha(i,j)-alpha_min)/alpha_min);
    end
    % cap & balance
    DE1(i,:) = max(min(DE1(i,:), max_shift*E0(i,:)), -max_shift*E0(i,:));
    DE1(i,:) = DE1(i,:) - mean(DE1(i,:));
end
L1     = sum(E0 - DE1,1);
profit1= sum(sum((repmat(rho,n,1)-alpha).*DE1));
DR1    = sum(DE1,1);

%% 3) Price-Driven Model (with Ep)
Ep = -0.3*eye(h) + 0.1*(ones(h)-eye(h));    % use your chosen Ep
pct_dev = (rho - base_price)./base_price;    % 1×h
DE2 = zeros(n,h);
for i = 1:n
    for j = 1:h
        DE2(i,j) = E0(i,j) * sum(Ep(j,:).*pct_dev);
    end
    DE2(i,:) = max(min(DE2(i,:), max_shift*E0(i,:)), -max_shift*E0(i,:));
    DE2(i,:) = DE2(i,:) - mean(DE2(i,:));
end
L2      = sum(E0 - DE2,1);
profit2 = sum(sum((repmat(rho,n,1)-alpha).*DE2));
DR2     = sum(DE2,1);

%% 4) Hybrid Model via combined elasticity
% Define Ea (incentive elasticity)
Ea =  0.8*eye(h) - 0.05*(ones(h)-eye(h));
% Combined elasticity
Ecomb = lambda*Ep + (1-lambda)*Ea;            % h×h
pct_dev = (rho - base_price)./base_price;     % 1×h
resp_vec = Ecomb * pct_dev.';                 % h×1

DE3 = zeros(n,h);
for i = 1:n
    for j = 1:h
        DE3(i,j) = E0(i,j) * resp_vec(j);
    end
    DE3(i,:) = max(min(DE3(i,:), max_shift*E0(i,:)), -max_shift*E0(i,:));
    DE3(i,:) = DE3(i,:) - mean(DE3(i,:));
end
L3      = sum(E0 - DE3,1);
profit3 = sum(sum((repmat(rho,n,1)-alpha).*DE3)) - gamma*sum(sum(DE3.^2));
DR3     = sum(DE3,1);

%% 5) Plot Load Profiles
figure; hold on;
plot(1:h,L_base,'k-','LineWidth',1.5);
plot(1:h,L1,    'b--','LineWidth',1.5);
plot(1:h,L2,    'g--','LineWidth',1.5);
plot(1:h,L3,    'r--','LineWidth',1.5);
xlabel('Hour','FontWeight', 'bold', 'FontSize', 12); ylabel('Load in kW','FontWeight', 'bold', 'FontSize', 12);
title('Demand Reduction Comparison','FontWeight', 'bold', 'FontSize', 12);
legend('Baseline','Incentive-NoElastic','Price-Elastic','Hybrid-Elastic','Location','Best','FontWeight', 'bold', 'FontSize', 12);
grid on;

%% 6) Profit Bar Chart with Values (Corrected)
figure;
profits = [profit1, profit2, profit3];  % Store profits in a vector
bh = bar(profits);
set(gca,'XTickLabel',{'Incentive Driven','Price Driven','Hybrid'}, 'FontSize',12,'FontWeight', 'bold');
ylabel('Profit (Rs)', 'FontSize',10);
title('ESP Profit Comparison', 'FontSize',12);
grid on;

% Calculate bar centers manually for MATLAB 2013
bar_centers = 1:length(profits); 

% Add text labels
for i = 1:3
    text(bar_centers(i), profits(i) + 3,...
        sprintf('Rs%.2f', profits(i)),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'FontSize',12,'FontWeight', 'bold');
end
%% 7) User-Level Profit (user 1)
u = 1;
up1 = sum((rho-alpha(u,:)).*DE1(u,:)) - gamma*sum(DE1(u,:).^2);
up2 = sum((rho-alpha(u,:)).*DE2(u,:)) - gamma*sum(DE2(u,:).^2);
up3 = sum((rho-alpha(u,:)).*DE3(u,:)) - gamma*sum(DE3(u,:).^2);
figure; bar([up1 up2 up3]);
set(gca,'XTickLabel',{'Incentive','Price','Hybrid'});
ylabel('Profit (Rs)');
title(sprintf('User %d Profit Contribution',u));
grid on;

%% 8) Baseline with Peak Marked (Modified)
figure;
plot(1:h,L_base,'k-o','LineWidth',2); 
hold on;

% Axis labels
xlabel('Hours', 'FontSize',12,'FontWeight', 'bold');
ylabel('Load in kW', 'FontSize',12,'FontWeight', 'bold');
title('Baseline Load with Peak Hours (19-22)', 'FontSize',14,'FontWeight', 'bold');
grid on;
set(gca,'XTick',1:h);
%% 9) Print Matrices & Stats
fprintf('\n--- PARAMETERS & MATRICES ---\n');
disp('E0 ='); disp(E0);
disp('rho ='); disp(rho);
disp('alpha ='); disp(alpha);
disp('Ep ='); disp(Ep);
disp('Ea ='); disp(Ea);
disp('Ecomb ='); disp(Ecomb);

fprintf('\n--- PROFIT REPORT ---\n');
fprintf('Incentive-NoElastic: Rs %.2f\n',profit1);
fprintf('Price-Elastic:       Rs %.2f\n',profit2);
fprintf('Hybrid-Elastic:      Rs %.2f\n',profit3);

fprintf('\n--- DEMAND REDUCTION STATS ---\n');
fprintf('Incentive max |?E|=%.2f at h=%s\n', max(DR1), mat2str(find(DR1==max(DR1))));
fprintf('Price     max |?E|=%.2f at h=%s\n', max(DR2), mat2str(find(DR2==max(DR2))));
fprintf('Hybrid    max |?E|=%.2f at h=%s\n', max(DR3), mat2str(find(DR3==max(DR3))));






