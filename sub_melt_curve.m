function [c1, c2] = sub_melt_curve(v, n, dv, norm_type)
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
%   v, n, dv = volume bins, counts, and bin spacing affected by sub. melt
%   norm_type = 2 for L2 and Inf for max norm
% OUTPUTS: FILL IN
% SYNTAX: FILL IN

% addpath('/Users/icebergs/general-code/') % add general functions to path
clf; % clear figure
season_cmap = 0.75*[0.000	0.000	0.000; 0.251	0.000	0.294; 0.463	0.165	0.514;
0.600	0.439	0.671; 0.761	0.647	0.812; 0.906	0.831	0.910;
0.851	0.941	0.827; 0.651	0.859	0.627; 0.353	0.682	0.380; 
0.106	0.471	0.216; 0.000	0.267	0.106; 0.50 0.50 0.50];

% call the fitting function
[error,c] = powerlaw_fit(v,n,norm_type)

% plot
loglog()


function [error, c] = powerlaw_fit(v,n,norm_type)
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
c = lsqcurvefit(@powerlaw_model, c0, v, n, [0 0], [1e12 1e12],opts); % fitting function with bounds
error = norm((powerlaw_model(c,v)-n)./n,norm_type);
end

end

