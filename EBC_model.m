function n = EBC_model(c,v) % Eq. 1 both sides
% Elastic Brittle fragmentation theory equation
% Eq. 1 from Astrom et al. 2021
% INPUTS:
%   c = vector of 5 model paramters
%   v = iceberg volume bins
% OUTPUTS:
%   n = modeled counts
% SYNTAX: n = EBC_model(c,v)

c1 = c(1); c2 = c(2); a = c(3); c3 = c(4); c4 = c(5); % parse input
n = c1 * v.^(-a) .* exp(-v/c2) + c3 * exp(-v/c4); % generate equation
end