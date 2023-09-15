function [OptionValue, stockTree, valueTree, deltaTree, gammaTree] = BinaryOptionTian(S0, X, r, sig, T, steps, earlyExercise)
% BinaryOptionTian: Price binary options using the Tian binomial model.
% Author: Samin Rajabi
% Syntax:
% [OptionValue, stockTree, valueTree, deltaTree, gammaTree] = ...
%    BinaryOptionTian(S0, X, r, sig, T, steps, earlyExercise)
%
% Inputs:
%   - S0: Initial stock price
%   - X: Strike price
%   - r: Risk-free interest rate
%   - sig: Volatility
%   - T: Time to maturity
%   - steps: Number of time steps to calculate
%   - earlyExercise: False for European, true for American
%
% Outputs:
%   - OptionValue: Calculated option value
%   - stockTree: Simulated potential future stock prices
%   - valueTree: Option prices
%   - deltaTree: Option sensitivity delta
%   - gammaTree: Option sensitivity gamma

% Validate input arguments
validateattributes(S0, {'numeric'}, {'positive', 'scalar'}, 'BinaryOptionTian', 'S0');
validateattributes(X, {'numeric'}, {'nonnegative', 'scalar'}, 'BinaryOptionTian', 'X');
validateattributes(r, {'numeric'}, {'scalar'}, 'BinaryOptionTian', 'r');
validateattributes(sig, {'numeric'}, {'positive', 'scalar'}, 'BinaryOptionTian', 'sig');
validateattributes(T, {'numeric'}, {'positive', 'scalar'}, 'BinaryOptionTian', 'T');
validateattributes(steps, {'numeric'}, {'positive', 'integer'}, 'BinaryOptionTian', 'steps');
validateattributes(earlyExercise, {'logical'}, {'scalar'}, 'BinaryOptionTian', 'earlyExercise');

% Calculate Tian binomial model parameters
dt = T / steps;
a = exp(r * dt);
V = exp(sig^2 * dt);
u = (0.5 * a * V) * (V + 1 + sqrt(V^2 + 2 * V - 3));
d = (0.5 * a * V) * (V + 1 - sqrt(V^2 + 2 * V - 3));
p = (a - d) / (u - d);

% Initialize stockTree
stockTree = zeros(steps + 1, steps + 1);
stockTree(1, 1) = S0;
for i = 2:steps + 1
    stockTree(1:i - 1, i) = stockTree(1:i - 1, i - 1) * u;
    stockTree(i, i) = stockTree(i - 1, i - 1) * d;
end

% Initialize valueTree
valueTree = zeros(size(stockTree));

% Calculate the value at expiry
inTheMoney = (stockTree(:, end) >= X);
valueTree(inTheMoney) = 1;
valueTree(~inTheMoney) = 0;

% Calculate option values
for j = steps:-1:1
    valueTree(1:j, j) = exp(-r * dt) * (p * valueTree(1:j, j + 1) + (1 - p) * valueTree(2:j + 1, j + 1));
    if earlyExercise
        intrinsic_value = double(stockTree(1:j, j) >= X);
        valueTree(1:j, j) = max(valueTree(1:j, j), intrinsic_value);
    end
end

% Calculate the option value
OptionValue = valueTree(1, 1);

% Calculate sensitivities
if nargout > 3
    dS = diff(stockTree);
    deltaTree = diff(valueTree) ./ dS;
    if nargout > 4
        gammaTree = diff(deltaTree) ./ diff(dS);
    end
end

