%% P219 - PS1 Q1.1

%% Question 1.1
clear
cd '/Users/willhotten/Library/CloudStorage/OneDrive-LondonBusinessSchool/Documents/PhD/Courses/Year 1/P219-Econometrics-I-Pt2/MATLAB'

%% Matrix preallocation

% Matrix to store average OLS estimates
mu_hat_mean_m = zeros(5,4);
phi_hat_mean_m = zeros(5,4);

% To reference position in loop over t
tref = 0;

%% Set loop for different values of phi and T
for t = [40, 80, 120, 280]

    % To reference position in loop over t and phi
    tref = tref + 1;
    phref = 0;
    
    for ph = [0.9, 0.95, 0.97, 0.99]

        % To reference position in loop over phi
        phref = phref + 1;

        %% Simulate AR(1)
        
        % Set parameter values
        mean_unc = 3;
        var_unc = 3;
        T = t;
        M = 5000;
        y0 = mean_unc;
        phi = ph;
        mu(1:T,1:M) = mean_unc*(1-phi);
        eps_mean = 0;
        eps_var = var_unc*(1-phi^2);
        
        % Set seed 
        rng(21)
        
        % Call AR1 function to simulate y
        [y] = AR1(T, M, y0, mu, eps_mean, eps_var, phi);
        
        %% Estimating AR(1)
        
        % Setting lag AR(p)
        p = 1;
        
        % Generating new X matrix with lagged y for OLS estimation
        X = [ones(T,1) mlag2(y,p)];
        
        % Creating OLS estimate matrices
        mu_hat = zeros(1,M);
        phi_hat = zeros(1,M);
        
        % OLS loop for each simulation
        for m = 1:M
            X_m = X(:,[1 m+1]);
            y_m = y(:,m);
            results = ols(y_m(p+1: end), X_m(p+1:end,:));
            mu_hat(1,m) = results.beta(1);
            phi_hat(1,m) = results.beta(2:end)';
        end
        
        % Calculate mean of OLS estimates
        mu_hat_mean = mean(mu_hat);
        phi_hat_mean = mean(phi_hat);

        % Store means for table
        mu_hat_mean_m(tref, phref) = mu_hat_mean;
        mu_hat_mean_m(5, phref) = mu(1,1);
        phi_hat_mean_m(tref, phref) = phi_hat_mean;
        phi_hat_mean_m(5, phref) = phi;
        
        
        %% Plotting OLS estimate distributions
        
        % Plots for mu
        path = '/Users/willhotten/Library/CloudStorage/OneDrive-LondonBusinessSchool/Documents/PhD/Courses/Year 1/P219-Econometrics-I-Pt2/MATLAB/Figures';

        % Distribution for mu plot
        mu_hat_tr = mu_hat';
        pd_mu = fitdist(mu_hat_tr,'kernel');
        %x_values_mu = -0.45:0.01:0.55;
        x_values_mu = 0:0.01:1;
        pdf_mu = pdf(pd_mu,x_values_mu);
        % Plotting
        figure(3), clf;
        area(x_values_mu,pdf_mu,'LineWidth',1.5,'EdgeColor','#0B0B45','FaceColor','#ADD8E6');
        xline(mu(t,m), 'Color', '#F47174', 'LineWidth', 1.75);
        xline(mu_hat_mean, 'Color', '#CC79A7', 'LineWidth', 1.75);
        grid on;
        hold on;
        temp=['Distribution of OLS estimates for mu, phi = ',num2str(phi),', T = ', num2str(T)];
        title(temp);
        legend('OLS estimates of mu','Actual value of mu', 'Average of OLS estimates');
        temp=[path,filesep,'dist_mu_phi_',num2str(phi),'_T_', num2str(T),'.png'];
        saveas(gca,temp);

        % Distribution for phi plot
        phi_hat_tr = phi_hat';
        pd_phi = fitdist(phi_hat_tr,'kernel');
        x_values_phi = 0.7:0.01:1.1;
        pdf_phi = pdf(pd_phi,x_values_phi);
        % Plotting
        figure(4), clf;
        area(x_values_phi,pdf_phi,'LineWidth',1.5,'EdgeColor','#0B0B45','FaceColor','#ACD1AF');
        xline(phi, 'Color', '#F47174', 'LineWidth', 1.75);
        xline(phi_hat_mean, 'Color', '#CC79A7', 'LineWidth', 1.75);
        grid on;
        hold on;
        temp=['Distribution of OLS estimates for phi, phi = ',num2str(phi),', T = ', num2str(T)];
        title(temp);
        legend('OLS estimates of phi','Actual value of phi', 'Average of OLS estimates');
        temp=[path,filesep,'dist_phi_phi_',num2str(phi),'_T_', num2str(T),'.png'];
        saveas(gca,temp);
    end
end

%% Create table of OLS estimates

% Specify column/row headers
VarNames = {'0.9', '0.95', '0.97', '0.99'};
RowNames = {'T = 40', 'T = 80', 'T = 120', 'T = 240', 'True'};

% Generate tables for mu and phi
tabmu = array2table(mu_hat_mean_m, 'RowNames',RowNames, 'VariableNames',VarNames);
tabphi = array2table(phi_hat_mean_m, 'RowNames',RowNames, 'VariableNames',VarNames);

% Export tables to LaTeX
temp=[path,filesep,'tab_mu'];
table2latex(tabmu, temp);
temp=[path,filesep,'tab_phi'];
table2latex(tabphi, temp);

