function [error,c1,c2] = powerlaw_fit(v, n, norm_type, lbounds, ubounds)
% This script fits a power law function to the iceberg
% size distribution data n(v) affected by submarine melt.
% The function is of the form:
% n(v) = c1*v^{-c2}
% and we are looking for parameters c = [c1 c2] given data n,v,dv.

% Author:  Jukes Liu
%           Department of Geosciences
%           Boise State University
%           Date: 04/07/2021
% INPUTS:
%   v, n = volume bins, counts affected by sub. melt
%   norm_type = 2 for L2, Inf for max norm, log for RMSLE
%   lbounds = vector of lower bounds on the 2 parameters
%   ubounds = vector of upper bounds on the 2 parameters
% OUTPUTS: FILL IN
% SYNTAX: FILL IN

addpath('/Users/jukesliu/Documents/CALVING-ICEBERGS/general/') % add general functions to path

% manually perform inversion to find initial guesses:
m = length(v);
A = [ones(m,1),v]; % initialize matrix A for linearization
pwr = log(n); % linearize n
b = pinv(A)*pwr; % Find model using Moore-Penrose Psuedoinverse
c0(1) = exp(b(1)); % c1
c0(2) = b(2); % c2 - power law exponent

c0 = abs(c0); %making sure initial guess is at least positive in all coeficients

% fit algorithm:
opts = optimset('Display','off');
c = lsqcurvefit(@powerlaw_model, c0, v, n, lbounds, ubounds ,opts); % fitting function with bounds

if norm_type == 'log' % if log norm is specified
    error = sqrt(mean((log10(powerlaw_model(c,v)+1) - log10(n+1)).^2));
else
    error = norm((powerlaw_model(c,v)-n)./n,norm_type);
end
c1 = c(1); c2 = c(2); % assign

end

