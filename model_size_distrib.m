%% model_size_distrib.m
% This script models iceberg size distributions stored in .mat structures
% using power laws and the Elastic-Brittle fragmentation theory curve,
% which corresponds to 'stable' iceberg calving
%
% Author:   Jukes Liu
%           Department of Geosciences
%           Boise State University
%           Date: 04/07/2021
function model_size_distrib(basepath,glacier_abbrev)
% clear all;
% addpath('/Users/ellynenderlin/mfiles/general/') % add general functions to path - change path

% set paths
% basepath='/Users/ellynenderlin/Research/NSF_Greenland-Calving/iceberg-calving/';
root_path = basepath; output_path = basepath;
cd_to_glacier = ['cd ''',root_path,'/',glacier_abbrev,'''']; 
eval(cd_to_glacier);
disp(['Analyzing ',glacier_abbrev]);

%% 0) grab the data
iceberg_dates = '';
% if length(dir([glacier_abbrev,'*_iceberg-data.mat'])) > 0
%     DEM_mats = dir([glacier_abbrev,'*_iceberg-data.mat']); % DEMs
% elseif length(dir([glacier_abbrev,'*_melange-DEMfilled.mat'])) > 0
    DEM_mats = dir([glacier_abbrev,'*_melange-DEMfilled.mat']); % DEMs
% end

for p = 1:numel(DEM_mats)
    fnames(p,:) = DEM_mats(p).name;
    iceberg_dates(p,:) = DEM_mats(p).name(4:11);
end
disp('Data identified');

nfiles = size(fnames,1); % grab number of files

nthresh = 1e-8; % set small size bin cutoff (n1 must be greater than this value)

% colormap:
season_cmap = 0.75*[0.000	0.000	0.000; 0.251	0.000	0.294; 0.463	0.165	0.514;
0.600	0.439	0.671; 0.761	0.647	0.812; 0.906	0.831	0.910;
0.851	0.941	0.827; 0.651	0.859	0.627; 0.353	0.682	0.380; 
0.106	0.471	0.216; 0.000	0.267	0.106; 0.50 0.50 0.50]; close gcf;

%% 1) plot size distributions to examine

for p = 1:nfiles
    % load the data
    name = fnames(p,:); 
    load_data = ['load ',name,' ''m''']; eval(load_data);
    disp(['Looping through iceberg distribution #',num2str(p),' of ',num2str(nfiles)]);
    v1 = double(m.melange.Asurfs(m.melange.bergs~=0)'); dv1 = double(m.melange.binwidth(m.melange.bergs~=0)'); 
    n1 = double(m.melange.bergs(m.melange.bergs~=0)')./dv1; % clear m;
    
    v1 = v1(n1>nthresh); n1 = n1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
    v1 = v1(~isnan(n1)); n1 = n1(~isnan(n1)); % Remove NaNs
    
    % plot
    figure(p); 
    loglog(v1,n1,'+','MarkerSize',6,'LineWidth',1,'Color',season_cmap(str2num(fnames(p,8:9)),:)); 
    grid on; set(gca,'FontSize',16); title(name(1:11)); %ylim([0.01,10^8]); 
    xlabel('surface area [m^2]'); ylabel('count');
    drawnow; disp('Plot generated.'); 
end

% %% 2A) Model-fitting with manual user input
% norm_type = 2;
% 
% for p = 1:nfiles
%     if p == 2;
%         % load the data
%         name = fnames(p,:);
%         load_data = ['load ',name]; eval(load_data);
%         disp(['Looping through iceberg distribution #',num2str(p),' of ',num2str(nfiles)]);
%         v1 = double(m.melange.Asurfs(m.melange.bergs~=0)'); dv1 = double(m.melange.binwidth(m.melange.bergs~=0)');
%         n1 = double(m.melange.bergs(m.melange.bergs~=0)')./dv1; % clear m;
%         
%         v1 = v1(n1>nthresh); n1 = n1(n1>nthresh); dv1 = dv1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
%         v1 = v1(~isnan(n1)); n1 = n1(~isnan(n1)); dv1 = dv1(~isnan(n1)); % Remove NaNs
%         
%         % plot
%         figure(p);
%         loglog(v1,n1,'+','MarkerSize',6,'LineWidth',1,'Color',season_cmap(str2num(fnames(p,8:9)),:)); hold on;
%         grid on; set(gca,'FontSize',16); title(name(1:11)); %ylim([0.01,10^8]);
%         xlabel('surface area [m^2]'); ylabel('count');
%         drawnow; disp('Plot generated.');
%         
%         % take user input
%         disp('Click to place the crossover point.');
%         a = ginput(1); 
%         x_cross = a(1); y_cross = a(2); % assign x,y from user input
%         disp('User input gathered.');
%         
%         % use input values and alpha = 2 as initial guesses
%         c2_guess = x_cross; % c2 = x coordinate from where curve deviates 
%         % from straight line in loglog space
%         c3_guess = y_cross; % c3 = y-intercept from the first term
%         a = 2; % hold alpha at 2
%         
%         % grab c1 and c4 guesses using EBC_fragmentation curve
%         [alpha,c1_guess,c2_x,c3_x,c4_guess,data_lims,error] = EBC_fragmentation_curve(name, v1, n1, dv1, norm_type, nthresh, 1, normalize_exp); % fit Eq. 1
%         
%         c0 = [c1_guess, c2_guess, a, c3_guess, c4_guess]; % initial guess vector
%         c0 = abs(c0); % make sure initial guess is at least positive in all coeficients   
%         disp(c0);
%         
%         % inversion with lsqcurvefit
%         opts = optimset('Display','off');
%         c = lsqcurvefit(@EBC_model, c0, v, n, [0 1e4 0 0 0], [1e12 1e12 3 1e12 1e12],opts); % ADJUST BOUNDS!!
%         disp(c);
%         
%         n_mod = EBC_model(c,v1); % grab model
%         
%         % add to plot
%         loglog(v1,n1,'+','MarkerSize',6,'LineWidth',1,'Color',season_cmap(str2num(name(8:9)),:)); hold on;
%         loglog(v1,n_mod,'LineWidth',2,'Color',season_cmap(str2num(name(8:9)),:)); hold on;
%         xlabel('surface area [m^2]'); ylabel('count'); set(gca,'FontSize',16);
%         title('Iceberg size distribution'); grid on;
%     end
% end

%% 2B) Automated model fitting
disp('Fitting the models...');
% INITIALIZE
nfiles = size(fnames,1); % grab number of files
for p = 1:nfiles; yrs(p) = str2num(fnames(p,4:7)); end % grab years
m = 0; bergs = [];
params = zeros(nfiles, 8); % matrix to hold filename and resulting parameters 
% c1, c2, c3, c4, alpha from E-BC, and c5, c6 from sub. melt power law 
% (will be NaNs if not used).

% MANUALLY SET THRESHOLDS:
use_dummy = 1; % if set to 1, the residuals calculated for determining 
% whether to account for submarine melt will be with respect to the dummy
% power law (approx. tangent but doesn't always work well). If not, the
% E-BC model will be used to calculate residuals.

taperthresh = 0.5; % the threshold value used to determine whether or not 
% counts taper due to submarine melt. Increase to increase the severity of 
% tapering required to be considered sub. melt influence.

dplawthresh = 10^5; % upper bound on the intercept for the dummy power 
% law fitting. Adjust to improve threshold.

norm_type = 2; % toggle between L2, max, and log norm using 2, Inf, and 'log'

normalize_exp = 1.2; % Increase to weight residuals towards end of the curve, minimum = 1

% LOOP START  
for p = 1:nfiles
%     if p > 21 % TOGGLE TO SELECT A SINGLE FILE FOR TESTING
        figure(p); set(gcf, 'Position', [500 500 1000 400]); % generate figure
        name = fnames(p,:); 
        load_data = ['load ',name,' ''m''']; eval(load_data);
        disp(['Looping through iceberg distribution #',num2str(p),' of ',num2str(nfiles)]);
        v1 = double(m.melange.Asurfs(m.melange.bergs~=0)'); dv1 = double(m.melange.binwidth(m.melange.bergs~=0)'); 
        n1 = double(m.melange.bergs(m.melange.bergs~=0)')./dv1; clear m;
        
        v1 = v1(n1>nthresh); n1 = n1(n1>nthresh); dv1 = dv1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
        v = v1(~isnan(n1)); n = n1(~isnan(n1)); % Remove NaNs
        
        % Fit E-BC fragmentation model to the data and plot
        [alpha,c1,c2,c3,c4,data_lims,error] = EBC_fragmentation_curve(name, v1, n1, dv1, norm_type, nthresh, 1, normalize_exp); % fit Eq. 1
        n_mod = EBC_model([c1,c2,alpha,c3,c4],v1); % grab model 

        subplot(1,2,1); % plot result
        loglog(v1,n1,'+','MarkerSize',6,'LineWidth',1,'Color',season_cmap(str2num(name(8:9)),:)); hold on;
        loglog(v1,n_mod,'LineWidth',2,'Color',season_cmap(str2num(name(8:9)),:)); hold on;
        xlabel('surface area [m^2]'); ylabel('count'); set(gca,'FontSize',16);
        title('Iceberg size distribution'); grid on;
        
        if use_dummy == 1
            % Fit dummy power law to determine if sub. melt is affecting the
            % smaller bergs in the size distribution:
            subset = 6:length(v1); % set subset of the distribution to fit 
            [error_d, c1_d, c2_d] = powerlaw_fit(v1(subset),n1(subset),norm_type,[max(n1) 0],[max(n1)*dplawthresh 1e12]); % fit with bounds, keep close to the first value
            nmod_dummy = powerlaw_model([c1_d, c2_d],v1);
            loglog(v1, nmod_dummy, '--', 'LineWidth',1,'Color',season_cmap(str2num(name(8:9)),:)); % plot the dummy power law
            res = n1 - nmod_dummy; % calculate residuals from dummy power law
        else % otherwise, use Eq. 1
            res = n1 - n_mod; % calculate residuals from E-BC model (data - model)
        end
            
        % Use negative residuals do determine if sub. melt is necessary
        resp = res(res >= 0); resn = res(res < 0); % split pos. and neg. residuals
        resp_idx = find(res >= 0); % grab index of positive residuals
        resn_idx = find(res < 0); % grab index of negative residuals
        
        subplot(1,2,2); % plot residuals
        loglog(resp_idx, resp./(n(resp_idx).^normalize_exp), 'rx'); hold on; 
        loglog(resn_idx, abs(resn)./(abs(n(resn_idx)).^normalize_exp), 'kx'); 
        legend('positive', 'negative'); set(gca,'FontSize',16);
        title('Residuals'); grid on;
        sgtitle(name(1:11)); % plot title
        
        submarine_melt_influence = 1; % initially assume influence of sub. melt
        % if there are at least 2 negative residuals in the first 8 points,
        % and the 1st or 2nd point has a negative residual:
        if (nnz(resn_idx < 8) > 3) && (resn_idx(1) <= 2)
            diffs = diff(resn_idx(resn_idx < 8)); % calculate differences sequentially
            [idxs] = find(diffs == 1); % find where the neg. residual indices are consecutive
            idx_sm_start = idxs(1)+1; % grab start index of sub. melt
            idx_sm_end = idxs(end)+2; % grab end index of sub. melt
            
            % EVALUATE WHETHER OR NOT THE DATA TAPER IN LOGLOG SPACE
            logres = log10(abs(res(idx_sm_start:idx_sm_end))); % log of the absolute value of residuals
            if nnz(diff(logres) < -1*taperthresh) > 0.5*length(logres) % if at least 50%
                % of the differences in log residal are above the threshold
                disp('Possible submarine melt influence discovered.');
            else; submarine_melt_influence = 0; end; % do not account for sub. melt 
        else; submarine_melt_influence = 0; end; % do not account for sub. melt
                  
        % Model submarine melt separately
        if submarine_melt_influence == 1 
            % separate data affected by sub. melt:
            v1_sm = v1(idx_sm_start:idx_sm_end);
            n1_sm = n1(idx_sm_start:idx_sm_end);
            dv1_sm = dv1(idx_sm_start:idx_sm_end);
            v1_EBC = [v1(1:idx_sm_start); v1(idx_sm_end:end)];
            n1_EBC = [n1(1:idx_sm_start); n1(idx_sm_end:end)];
            dv1_EBC = [dv1(1:idx_sm_start); dv1(idx_sm_end:end)];
            
            % fit power law model to data influenced by submarine melt
            [error2,c5,c6] = powerlaw_fit(v1_sm, n1_sm, 2,[0 0],[1e12 1e12]);
            nmod_sm = powerlaw_model([c5,c6],v1(1:8));
            
            if c1 < c5 % if the sub. melt y-intercept is greater than the original model y-intercept
                disp('No submarine melt.');
            else % otherwise,
                disp('Submarine melt influence confirmed.')
                
                % fit E-BC fragmentation curve to data not influenced by sub. melt
                [alpha,c1,c2,c3,c4,data_lims,error] = EBC_fragmentation_curve(name, v1_EBC, n1_EBC, dv1_EBC, norm_type, nthresh, 1, normalize_exp); % fit Eq. 1
                nmod_EBC = EBC_model([c1,c2,alpha,c3,c4],v1_EBC);

                % plot over 2nd subplot
                subplot(1,2,1); % old result
                loglog(v1,n1,'+','MarkerSize',6,'LineWidth',1,'Color',season_cmap(str2num(name(8:9)),:)); hold on;
                loglog(v1,n_mod,'LineWidth',2,'Color',season_cmap(str2num(name(8:9)),:)); hold on; % original model
                if use_dummy == 1
                    loglog(v1, nmod_dummy, '--', 'LineWidth',2,'Color',season_cmap(str2num(name(8:9)),:)); %  dummy power law
                end
                xlabel('surface area [m^2]'); ylabel('count'); set(gca,'FontSize',16);
                title('Iceberg size distribution'); grid on;
                subplot(1,2,2); % new result
                loglog(v1,n1,'+','MarkerSize',6,'LineWidth',1,'Color',season_cmap(str2num(name(8:9)),:)); hold on;
                loglog(v1(1:8), nmod_sm,'k--','LineWidth',1); hold on; grid on;
                loglog(v1_EBC,nmod_EBC,'LineWidth',2,'Color',season_cmap(str2num(name(8:9)),:)); hold on;
                xlabel('surface area [m^2]'); ylabel('count'); set(gca,'FontSize',16);
                title('submarine melt considered');
                legend('data','submarine melt','E-BC frag. model');
                sgtitle(name(1:11)); % plot title
            end
        else
            c5 = NaN; c6 = NaN;
            disp('Distribution not influenced by submarine melt.')
        end
        
        % save the figure
        if not(isfolder([basepath,'/',glacier_abbrev,'/models/']))
            mkdir([basepath,'/',glacier_abbrev,'/models/']) % make models folder if it doesn't exist
            disp('Models folder created.');
        end
        saveas(gcf,[basepath,glacier_abbrev,'/models/',name(1:11),'_model.png']) % save into models folder
        
        % save the parameters
        params(p,:) = [str2num(name(4:11)), c1, c2, alpha, c3, c4, c5, c6];
        
        % clear variables
        clearvars -except nfiles fnames season_cmap nthresh basepath root_path glacier_abbrev output_path use_dummy taperthresh dplawthresh params norm_type normalize_exp
        
%     end % TOGGLE TO SELECT A SINGLE FILE FOR TESTING
end

% export parameters to a csv file
dt = string(datetime('now')); dt = strrep(dt,' ','_');
dt = strrep(dt,':',''); dtstring = extractBefore(dt(1),17); % grab datetime
delete_csvs = ['delete ',strcat(basepath,glacier_abbrev,'/models/',glacier_abbrev,'_parameters_*.csv')]; eval(delete_csvs);
writematrix(params,strcat(basepath,glacier_abbrev,'/models/',glacier_abbrev,'_parameters_',dtstring,'.csv'));

%% FUNCTIONS

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