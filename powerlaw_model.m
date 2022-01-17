function n = powerlaw_model(c,v)
% Submarine melting power law
% INPUTS:
%   c = vector of 2 model paramters
%   v = iceberg volume bins
% OUTPUTS:`
%   n = modeled counts
% SYNTAX: n = powerlawfit(c,v)

c1 = c(1); c2 = c(2); % parse input
n = c1 * v.^(-c2); % generate equation
end