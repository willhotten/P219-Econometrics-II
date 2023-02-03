%% P219 - PS1 Q1.2

%% Question 1.2
clear
cd '/Users/willhotten/Library/CloudStorage/OneDrive-LondonBusinessSchool/Documents/PhD/Courses/Year 1/P219-Econometrics-I-Pt2/MATLAB'


%% Setup

% Set parameter values
mean_unc(1:100,:) = 5;
mean_unc(101:192,:) = 0;
mean_unc(193:228,:) = 4.5;
mean_unc(229:280,:) = 0.5;
var_unc = 3;
T = 280;
M = 1;
y0 = mean_unc(1);
phi = 0.9;
mu = mean_unc*(1-phi);
eps_mean = 0;
eps_var = var_unc*(1-phi^2);

% Set seed (8 specified in assignment)
rng(8)

%% Generate samples series for y to be called ytrue
[ytrue] = AR1(T, M, y0, mu, eps_mean, eps_var, phi);

% Plot the generated series and unconditional mean
path = '/Users/willhotten/Library/CloudStorage/OneDrive-LondonBusinessSchool/Documents/PhD/Courses/Year 1/P219-Econometrics-I-Pt2/MATLAB/Figures';

figure(1), clf;
plot(ytrue, 'b');
hold on;
plot(mean_unc, 'Color', 'r');
hold off;
grid on;
%set(gca,'Yticklabel',[]);
xlabel('Quarters from 1948Q1');
ylabel('y');
%temp=['Plot of ytrue and E[y_t]'];
%title(temp);
legend('y true', 'E[y_t]');
temp=[path,filesep,'ts_y_true','.png'];
saveas(gca,temp);

%% Setting parameters and preparing for various forecasts

% Setting lag AR(p), window length (h), and end date of data examined (s)
p = 1;
h = 12;
s = 168; % Corresponds to 1989Q4

% Generating X matrix with lagged y for AR(p) estimation
X = [ones(T,1) mlag2(ytrue,p)];

% Create symbolic variable for summations that follow
syms j

%% Estimating expanding window - researcher 1

% Pre-allocate forecast matrix
yf_1 = zeros(h,T-s);

% Loop for expanding window
for w = 0:(T-s-1)

    % Extracting subset window of ytrue and placing in matrix
    ytrue_window = ytrue(1:s+w);
    X_window = X((1:s+w),:);
    
    % OLS computation
    results = ols(ytrue_window(p+1: end), X_window(p+1:end,:));
    mu_hat = results.beta(1);
    phi_hat = results.beta(2:end)';
    
    % Loop to generate 12 period ahead forecast
    for i=1:h
        yf_1(i,w+1) = mu_hat*symsum(phi_hat^(j),j,0,i) + (phi_hat^i)*ytrue_window(end);
    end

end

%% Estimating rolling window - researcher 2

% Set lag size of rolling window (= no. of quarters - 1)
size = 39;

% Pre-allocate forecast matrix
yf_2 = zeros(h,T-s);

% Loop for rolling window
for w = 0:(T-s-1)

    % Extracting subset window of ytrue and placing in matrix
    ytrue_window = ytrue(s+w-size:s+w);
    X_window = X((s+w-size:s+w),:);
    
    % OLS computation
    results = ols(ytrue_window(p+1: end), X_window(p+1:end,:));
    mu_hat = results.beta(1);
    phi_hat = results.beta(2:end)';
    
    % Loop to generate 12 period ahead forecast
    for i=1:h
        yf_2(i,w+1) = mu_hat*symsum(phi_hat^(j),j,0,i) + (phi_hat^i)*ytrue_window(end);
    end

end

%% Random walk forecast - researcher 3
% Forecasting a random walk: E[yt+1] = yt and so on

% Pre-allocate forecast matrix
yf_3 = zeros(h,T-s);

% Forecast loop
for w = 0:(T-s-1)

    % Select final observation available of y for forecast
    yfinal = ytrue(s+w);
    
    % Loop to generate 12 period ahead forecast
    for i=1:h
        yf_3(i,w+1) = yfinal;
    end

end

%% Forecasts using true model - researcher 4

% Pre-allocate forecast matrix
yf_4 = zeros(h,T-s);

% Forecast loop
for w = 0:(T-s-1)

    % Select final observations of y and mu
    yfinal = ytrue(s+w);
    mu_final = mu(s+w);
    
    % Store values of mu for time periods in forecast window
    mu_window = zeros(h,1);
    for i=1:h
        if s+w+i <= T
                mu_window(i) = mu(s+w+i);
            else
                mu_window(i) = mu_final;
        end
    end
      
    % Generating forecast knowing true parameter values
    for i = 2:h
        yf_4(1,w+1) = mu_window(1) + phi*yfinal;
        yf_4(i,w+1) = mu_window(i) + phi*yf_4(i-1,w+1);
    end

end

%% Constructing comparable series for ytrue

% Pre-allocate yf_true matrix
yf_true = zeros(h,T-s);

% Forecast loop
for w = 0:(T-s-1)

% Loop for yf_true
    for i=1:h      
        if s+w+i <= T
            yf_true(i,w+1) = ytrue(s+w+i);
        else
            yf_true(i,w+1) = NaN;
        end
    end

end

%% Comparing forecasts to ytrue

% Preallocate RMSE and MAE
RMSE = zeros(4,h);
MAE = zeros(4,h);

% Loop to calculate RMSE and MAE for each forecaster and window (can be made more efficient)
for w = 1:h
    
    RMSE(1,w) = sqrt(nanmean((yf_true(w,:) - yf_1(w,:)).^2));
    MAE(1,w) = nanmean(abs(yf_true(w,:) - yf_1(w,:)));
    RMSE(2,w) = sqrt(nanmean((yf_true(w,:) - yf_2(w,:)).^2));
    MAE(2,w) = nanmean(abs(yf_true(w,:) - yf_2(w,:)));
    RMSE(3,w) = sqrt(nanmean((yf_true(w,:) - yf_3(w,:)).^2));
    MAE(3,w) = nanmean(abs(yf_true(w,:) - yf_3(w,:)));
    RMSE(4,w) = sqrt(nanmean((yf_true(w,:) - yf_4(w,:)).^2));
    MAE(4,w) = nanmean(abs(yf_true(w,:) - yf_4(w,:)));

end

%% Plotting

% Plotting RMSE ratios
figure(2), clf;
plot(RMSE(1,:)./RMSE(4,:),'LineWidth',1);
hold on;
plot(RMSE(2,:)./RMSE(4,:),'LineWidth',1);
hold on;
plot(RMSE(3,:)./RMSE(4,:),'LineWidth',1);
hold off;
grid on;
xlabel('Window length');
ylabel('RMSE ratio');
%temp=['Plot of RMSE'];
%title(temp);
legend('Researcher 1', 'Researcher 2', 'Researcher 3');
temp=[path,filesep,'RMSE','.png'];
saveas(gca,temp);

% Plotting MAE ratios
figure(3), clf;
plot(MAE(1,:)./MAE(4,:),'LineWidth',1);
hold on;
plot(MAE(2,:)./MAE(4,:),'LineWidth',1);
hold on;
plot(MAE(3,:)./MAE(4,:),'LineWidth',1);
hold off;
grid on;
xlabel('Window length');
ylabel('MAE ratio');
%temp=['Plot of MAE'];
%title(temp);
legend('Researcher 1', 'Researcher 2', 'Researcher 3');
temp=[path,filesep,'MAE','.png'];
saveas(gca,temp);