

% Author:Samin Rajabi 
% Function to price an American Put option using the Crank-Nicolson Scheme.
clc 
clear all

function [price, opt_cost, exitflag] = AmericanPut_CrankNicolson(Svec, X, r, sig, tvec)
    N = length(tvec) - 1;
    M = length(Svec) - 1;
    dt = tvec(2) - tvec(1);

    % Initialize matrices
    price = zeros(M+1, N+1);
    intrinsic_value = max(X - Svec, 0);
    price(:, N+1) = intrinsic_value;
    price(1, :) = X;
    price(M+1, :) = 0;

    alpha = (dt/4) * (sig^2 * (1:M).^2 - r * (1:M));
    beta = -(dt/2) * (sig^2 * (1:M).^2 + r);
    gamma = (dt/4) * (sig^2 * (1:M).^2 + r * (1:M));

    C = diag(-alpha(3:M), -1) + diag(1 - beta(2:M)) + diag(-gamma(2:M-1), 1);
    D = diag(alpha(3:M), -1) + diag(1 + beta(2:M)) + diag(gamma(2:M-1), 1);

    opt_cost = nan(1, N);
    exitflag = nan(1, N);

    opts = optimoptions('lsqnonlin');
    opts.Display = 'none';
    opts.MaxIterations = 5000;
    opts.FunctionTolerance = 1e-12;
    opts.OptimalityTolerance = 1e-12;
    opts.StepTolerance = 1e-12;

    for idx = N:-1:1
        offset = [alpha(2) * (price(1, idx) + price(1, idx+1));
                  zeros(M-2, 1);
                  gamma(M-1) * (price(M+1, idx) + price(M+1, idx+1))];

        F0 = price(2:M, idx+1);
        rhs = D * F0 + offset;

        costFcn = @(F) (C * F - rhs);
        [Fopt, opt_cost(idx), ~, exitflag(idx)] = lsqnonlin(costFcn, F0, [], [], opts);

        price(2:M, idx) = Fopt;
    end
end
