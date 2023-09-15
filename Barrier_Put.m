% Script to price a barrier option using a monte-carlo approach
% Author: Samin Rajabi
% Created: 2019

clear all
close all
clc

S0 =50;
Sb = 48;
X = 55;
r = 0.04;
sig = 0.1;
T = 1;

nruns = 10000;

% This is not a path dependent option so only the final value is required
tic
nu = r - sig*sig/2;
Sf = S0*exp(nu*T+sig*sqrt(T)*randn(1,nruns));
% factor in the put
payoffT = max(X-Sf,0);
% Include the barrier
payoffT(payoffT>(X-Sb))=0;

str = ['Monte-Carlo pricing a simple barrier/put option with (',num2str(nruns),' paths).'];
disp(str);
mean(payoffT)*exp(-r*T)
toc
