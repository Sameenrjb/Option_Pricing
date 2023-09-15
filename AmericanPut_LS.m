% Script to price an American PUT using the least squares approach of
% Longstaff and Schwartz, 
% Author: Samin Rajabi
% a continuation boundary is calculated using one set of simulations and
% then the option price is calculated using a second set of simulations.

clear all
close all

rng('default');

% Define parameters
S0 = 50;
X = 55;
sig = 0.2;
r = 0.04;
dt = 1/365;
steps = 365;
nPaths = 10000;
oType = 'put';

% Since we're plotting a potentially very large number of paths, create a
% figure that is invisible, then if we want to see the plot we can just
% make it visible
hf = figure('Visible','off');

% Generate Asset Paths
S = AssetPaths(S0,r,sig,dt,steps,nPaths);

% Calculate the discount factor for each time step
disc = exp(-r*dt);

% Calculate intrinsic values (these correspond to potential cashflows)
cFlows = max(X-S,0);
% Preallocate an array for the coefficients
cMat = nan(steps,3);

% Loop backwards in time filling the cash flows matrix
for idx = steps:-1:2
    % Determine if there are any positive cashflows to regress
    mask = (cFlows(idx,:) > 0);
    % ehich cash flows are greater than zero
    if any(mask)
        % form the Y and X columns to be regressed
        Xdata = S(idx,mask);
        Ydata = (disc.^(1:(steps-idx+1)))*cFlows(idx+1:end,mask);
        % Do the regression, in this case using a quadratic fit (if
        % possible)
        if numel(Xdata) == 1
            coeffs = Ydata;
        elseif var(Xdata) < eps
            coeffs = 0;  % There is no variation to regress
        else
            coeffs = polyfit(Xdata,Ydata,min(2,numel(Xdata)-1));
        end
        % Calculate potential payoff from "backwards induction"
        payoffs = polyval(coeffs,Xdata);
        % Find location(s) where early exercise is optimal
        eeLoc = (cFlows(idx,mask) > payoffs);
        % Update the cash flow matrix to account for early exercise
        mask(mask) = eeLoc;  % These locations get exercised early
        cFlows(idx,:) = mask.*cFlows(idx,:);  % These are the cashflows
        cFlows((idx+1):end,mask) = 0;  % zero everything after the early ex.
        % Save the basis functions (for use in the extended method)
        cMat(idx,:) = coeffs;
    end
end

% Generate a new set of simulations and plot them
S = AssetPaths(S0,r,sig,dt,steps,nPaths);
% generate a plot which is used to cycle through the time points and show
% the continuation boundary at each time
payoff = 0; % value to accumulate the discounted payoffs
for idx = 2:steps
    coeffs = cMat(idx,:);
    % Plot the simulated points at this time
    simX = S(idx,:);
    simY = max(X-simX,0);  % Intrinsic value
    itmLoc = (simX<X);  % For a put
    contVal = zeros(size(simY));
    contVal(itmLoc) = polyval(coeffs,simX(itmLoc)); % Continuation value
    contLoc = (simY > contVal);
    if any(contLoc)
        % Set the rest of this path to Nan so it doesn't plot
        S(idx+1:end,contLoc) = nan;
        % Plot a point at the end of the line to highlight it's end
        line(repmat(idx,sum(contLoc),1),simX(contLoc),...
            'Marker','.','LineStyle','none',...
            'MarkerSize',20,'Color','r');
        hold on
        % Calculate/accumulate the discounted payoff for this path
        payoff = payoff + exp(-r*dt*idx)*sum(max(X-S(idx,contLoc),0));
    end
end
% Plot the paths
plot(S);
title(sprintf('%d Simulation Paths with Continuation Boundary',nPaths));
grid on;
line([0 steps+1],X*[1 1],'Color','r','LineWidth',2);

% Calculate/accumulate the discounted payoff for paths that made it to
% expiry
payoff = payoff + exp(-r*dt*steps)*sum(max(X-S(end,~isnan(S(end,:))),0));
oPrice = payoff/nPaths;

% Display the price
fprintf('The calculated price using the ''extended'' approach is: %.4f\n',...
    oPrice);

% Make the figure large and then make it visible
fullsize = get(0,'ScreenSize');
set(hf,'OuterPosition',[0 1 0.9*fullsize(3:4)]);
movegui(hf,'center');
set(hf,'Visible','on');

