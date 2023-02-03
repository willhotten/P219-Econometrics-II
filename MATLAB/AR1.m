%%% AR(1) function

%% Function setup
function [y] = AR1(T, M, y0, mu, eps_mean, eps_var, phi);

% Preallocate y and epsilon matrices
y = zeros(T,M);
eps = zeros(T,M);

% Simulate AR(1)
for m = 1:M

    for t = 2:T
        % Set initial value for y using specified y0
        y(1,m) = y0;
    
        % Simulate value for epsilon from normal distribution of epsilon
        eps(t,m) = normrnd(eps_mean, sqrt(eps_var));
    
        % Iterate value for y using AR1 formula
        y(t,m) = phi*y(t-1,m) + mu(t,m) + eps(t,m);

    end
    
end

