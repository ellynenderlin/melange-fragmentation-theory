function [alpha,c1,c2,c3,c4,data_lims,error] = EBC_fragmentation_curve(fname, v1, n1, dv1,norm_type, nthresh, tracktime, normalize_exp)
% This script fits the fragmentation theory function to the iceberg
% size distribution data n(v). The function is of the form:
% n(v) = c1*v^{-a}*exp(-v/c2) + c3*dv*exp(-v/c4)
% and we are looking for parameters c = [c1 c2 a c3 c4]'
% given data n,v,dv.
%
% Authors:  Jukes Liu
%           Department of Geosciences
%           Boise State University
%           Date: 04/07/2021
%
%           Michal A. Kopera
%           Department of Mathematics
%           Boise State University
%           Date: 05/12/2020

% INPUTS:
%   fname = filename
%   v1, n1, dv1 = volume bins, counts, and bin spacing for the size distribution
%   norm_type = specify the norm: 2 = L2 norm, Inf = maximum norm, log = RMSLE
%   nthresh = threshold on count to remove small fragments
%   tracktime = 0 to not display processing time, 1 to display
%   normalize_exp = exponent to modulate normalization of data
% OUTPUTS: FILL IN
% SYNTAX: FILL IN

%initialize
addpath('/Users/jukesliu/Documents/CALVING-ICEBERGS/general/') % add general functions to path
clf; 
if tracktime == 1; tic; end % track time
season_cmap = 0.75*[0.000	0.000	0.000; 0.251	0.000	0.294; 0.463	0.165	0.514;
0.600	0.439	0.671; 0.761	0.647	0.812; 0.906	0.831	0.910;
0.851	0.941	0.827; 0.651	0.859	0.627; 0.353	0.682	0.380; 
0.106	0.471	0.216; 0.000	0.267	0.106; 0.50 0.50 0.50];

if sum((n1.*dv1).*v1) >= 175e3 % if DEM coverage is substantial
    
    %remove the very small iceberg fragments of BIG icebergs since those
    %are simply elevation peaks in the largest icebergs
    v1 = v1(n1>nthresh); n1 = n1(n1>nthresh);
    
    v = v1(~isnan(n1)); n = n1(~isnan(n1)); % NO REMOVAL OF FIRST POINTS
    
    %begin the search loop - compute error for every combination of i_cut, j_cut
    err_array = 1e12*ones(length(v)); %pre-allocate error array
    c_array = zeros(length(v),length(v),5);
%     a_array = zeros(length(v),length(v));
    
    for i=4:length(v)-3
        for j=4:length(v)-3
            %create a function which computes error for parameter a
            %given the cutoff values i_cut and j_cut
            %fun = @(a) fit(v,n,i_cut(i),j_cut(j),a,norm_type);
            
            %find for which a the error is minimal
            %a = fminbnd(fun,1.5, 2.0);
            a = 2;
            %find parameters and error for the fit with new a
            [error, c] = EBC_fit(v,n,i,j,a,norm_type);
            
            err_array(i,j) = error;
            c_array(i,j,:) = c;
%             a_array(i,j) = a;
        end
    end
    
    %find global minimum in the i_cut, j_cut landscape
    %basically looking for which cutoff values we get minimal error
    [min_val,idx]=min(err_array(:));
    
    [i,j]=ind2sub(size(err_array),idx);
    c = c_array(i,j,:);
%     a = a_array(i,j);
    
    %rename variables
    alpha=a;
    c1 = c(1); c2 = c(2); c3 = c(4); c4 = c(5);
    data_lims(:) = [i j]; error = min_val;
    
    %display results
    fprintf(fname(4:11)); fprintf("\n");
    fprintf("[c1, c2, c3, c4] = [%e %e %e %e]\n",c(1),c(2),c(4),c(5));
    fprintf("alpha = %f\n",c(3));
    fprintf("==========================================\n")
      
    if tracktime == 1; toc; end
    
else
    fprintf("==========================================\n")
    fprintf(fname(4:11)); fprintf("\n");
    fprintf("DEM covers too small of an area to get reliable data\n")
    fprintf("==========================================\n")
    
    %fill-in NaNs as outputs
    alpha=NaN;
    c1 = NaN; c2 = NaN; c3 = NaN; c4 = NaN;
    data_lims(:) = [NaN NaN]; error = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit function:
function [error,c] = EBC_fit(v,n,i_cut,j_cut,a,norm_type)
%the actual fitting function

%exponential fit - left side
vl = v(1:i_cut); nl = n(1:i_cut);

m = length(vl);
A = [ones(m,1),vl];
lhs = log(nl) - log(vl.^(-a)); % linearize by taking hte log
% b = (A'*A)\(A'*lhs);
b = pinv(A)*lhs; % psuedoinverse
cl(1) = exp(b(1));
cl(2) = -1/b(2);
cl(3) = a;

%exponential fit - right side
vr = v(j_cut:end); nr = n(j_cut:end);

m = length(vr);
A = [ones(m,1),vr];
rhs = log(nr); % linearize by taking the log
% b = (A'*A)\(A'*rhs);
b = pinv(A)*rhs; % psuedoinverse
cr(1) = exp(b(1));
cr(2) = -1/b(2);

%now use the linear approach as an initial guess to fit the entire function
c0 = [cl, cr]; %Initial guess
c0 = abs(c0); %making sure initial guess is at least positive in all coeficients

opts = optimset('Display','off');
c = lsqcurvefit(@EBC_model, c0, v, n, [0 1e4 0 0 0], [1e12 1e12 3 1e12 1e12],opts); % ADJUST BOUNDS!!

if norm_type == 'log' % if log norm is specified
    error = sqrt(mean(((log10(EBC_model(c,v)+1) - log10(n+1))./(n.^normalize_exp)).^2)); % RMSLE
else % otherwise use L2 or Inf in the norm function
    error = norm((EBC_model(c,v)-n)./(n.^normalize_exp),norm_type); % normalized residuals
end
end

end