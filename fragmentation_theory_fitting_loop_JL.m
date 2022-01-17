function [alpha,c1,c2,c3,c4,c5,beta,data_lims,errors] = fragmentation_theory_fitting_loop_JL(fnames,iceberg_dates)
%% This script fits the fragmentation theory function to the iceberg
% size distribution data n(v). The function is of the form:
% n(v) = c1*v^{-a}*exp(-v/c2) + c3*dv*exp(-v/c4)
% and we are looking for parameters c = [c1 c2 a c3 c4]'
% given data n,v,dv.
% The script uses a brute force approach to searching, which could be
% improved with steepest descent search instead
%
% Author:   Michal A. Kopera
%           Department of Mathematics
%           Boise State University
% Date:     05/12/2020
% Modified by:  Jukes Liu
%               Department of Geosciences
% Date:         03/24/2021

%% For testing:
% set paths
clear all;
basepath='/Users/jukesliu/Documents/CALVING-ICEBERGS/';
glacier_abbrev = 'KO';
root_path = basepath; output_path = basepath;
addpath('/Users/jukesliu/Documents/CALVING-ICEBERGS/general/') % add general functions to path

cd_to_folder = ['cd ',basepath,'/',glacier_abbrev,'/']; eval(cd_to_folder); % cd to the folder
DEM_mats = dir([glacier_abbrev,'*_melange-DEMfilled.mat']); iceberg_dates = ''; %DEMs
% DEM_mats = dir([glacier_abbrev,'*_iceberg-data.mat']); iceberg_dates =
% ''; % iceberg data
% for i = 1:length(DEM_mats)
%     iceberg_dates(i,:) = DEM_mats(i).name(4:11);
% end
for p = 1:numel(DEM_mats)
    fnames(p,:) = DEM_mats(p).name;
end

%% Actual code:

% 1) INITIALIZE
clf; %tic

%choose error norm to use:
% 2 - L2 norm (default)
% Inf - maximum norm
norm_type = 2;
% 
% %choose bin cutoff size to account for submarine melt-induced deviation
% %from a power law for the smallest icebergs
% berg_cutoff = 1;

%prepare a figure
% fitplot = figure; 
for p = 1:size(fnames,1); yrs(p) = str2num(fnames(p,4:7)); end
% set(fitplot,'position',[50 50 500 ceil(length(unique(yrs))/2)*300]);
%pull indices for each year
% time_cmap = colormap(parula(length(fnames)));
season_cmap = 0.75*[0.000	0.000	0.000; 0.251	0.000	0.294; 0.463	0.165	0.514;
0.600	0.439	0.671; 0.761	0.647	0.812; 0.906	0.831	0.910;
0.851	0.941	0.827; 0.651	0.859	0.627; 0.353	0.682	0.380; 
0.106	0.471	0.216; 0.000	0.267	0.106; 0.50 0.50 0.50]; close gcf;

% read data
pl = []; date_names = []; ref = 0;
nfiles = size(fnames,1);

% SUBSET DATA FOR TESTING

%% Start loop
for p = 1:nfiles % loop through files
    
    if fnames(p,:) == 'KO-20150419_melange-DEMfilled.mat' % specify a file for testing
        
        load_data = ['load ',fnames(p,:)]; eval(load_data); % load the data
        disp(['Looping through iceberg distribution #',num2str(p),' of ',num2str(nfiles)]);

        % grab iceberg size distribution and counts
        v1 = double(m.melange.Asurfs(m.melange.bergs~=0)'); dv1 = double(m.melange.binwidth(m.melange.bergs~=0)'); 
        n1 = double(m.melange.bergs(m.melange.bergs~=0)')./dv1; clear m;

        % establish plot variables
        if p == find(yrs == str2num(fnames(p,4:7)),1,'first'); ref = 0; date_names = []; figure; end
        ref = ref+1; date_names = [date_names; fnames(p,4:11)];

        % if DEM is large enough to get accurate coverage:
        if sum((n1.*dv1).*v1) >= 175e6

            %remove the very small iceberg fragments of BIG icebergs since those
            %are simply elevation peaks in the largest icebergs
            v1 = v1(n1>0.01); n1 = n1(n1>0.01);

            %remove the first two points (influenced by submarine melting)
            %we may need a more automatic way of choosing how many points to
            %neglect
    %         v = v1(berg_cutoff:end); n = n1(berg_cutoff:end); % no cutoff
    %         v = v(~isnan(n)); n = n(~isnan(n));

            v = v1(~isnan(n1)); n = n1(~isnan(n1)); % NO REMOVAL OF FIRST POINTS

            %begin the search loop - compute error for every combination of
            %i_cut, j_cut, k_cut
            err_array = 1e12*ones(length(v),length(v),length(v)); %pre-allocate error array with high error values 
            % that would not influence minsearch
            c_array = NaN*zeros(length(v),length(v),length(v), 7); % 3 variables to search through, 7 parameters c
            a_array = NaN*zeros(length(v),length(v),length(v)); % 3 variables to search through

            for k=1:round(length(v)/2) % search the first half of the dataset
                for i=1:length(v)-3
                    for j=1:length(v)-3           
                        % Find for which a the error is minimal
                        a = 2; % set constant with value based on fragmentation theory

                        % Find parameters and error for the fit with new a
                        [error, c] = fit(v,n,k,i,j,a,norm_type);

                        err_array(i,j,k) = error;
                        c_array(i,j,k,:) = c(1:7); % don't keep c(8), which is just k
                        a_array(i,j,k) = a;
                    end
                end
            end

            %find global minimum in the i_cut, j_cut landscape
            %basically looking for which cutoff values we get minimal error
            [min_val,idx]=min(err_array(:)); % find minimum error
            [i,j,k]=ind2sub(size(err_array),idx); % grab the i, j, k associated with it
            c = c_array(i,j,k,:); % grab the parameters
            a = a_array(i,j,k); % grab the exponent a, should be the constant set earlier

            % Rename variables to match the equations
            alpha(p)=c(3);
            c1(p) = c(1); c2(p) = c(2); c3(p) = c(4); c4(p) = c(5); 
            c5 = c(6); c6(p) = c(7);
            data_lims(p,:) = [i j k]; errors(p) = min_val;

            % Display results
            fprintf(fnames(p,4:11)); fprintf("\n");
            fprintf("[c1, c2, c3, c4, c5] = [%e %e %e %e %e]\n",c1,c2,c3,c4,c5);
            fprintf("alpha = %f\n",alpha);
            fprintf("c6 = %f\n",c6);
            fprintf("==========================================\n")
            %     toc

            c = [c1, c2, alpha, c3, c4, c5, c6, k]; % arrange c and add k as a parameter for plotting
            % Plot result
            pl(ref) = loglog(v1,n1,'s','MarkerSize',10,'LineWidth',1,'Color',season_cmap(str2num(fnames(p,8:9)),:)); hold on;
            loglog(v1,model(c,v1),'LineWidth',2,'Color',season_cmap(str2num(fnames(p,8:9)),:)); hold on;
            xlabel('surface area [m^2]'); ylabel('count');
            ylim([0.01,10^8]); %ylim([0.01,10*max(n1)]);
    %         if p == find(yrs == str2num(fnames(p,4:7)),1,'last')
    %             legend(pl,date_names); 
    %             if isinf(norm_type)
    %                 saveas(gcf,[fnames(p,1:2),'-iceberg-size-distributions_',fnames(p,4:7),'-',num2str(berg_cutoff),'mcutoff-INFnorm.png'],'png'); 
    %             else
    %                 saveas(gcf,[fnames(p,1:2),'-iceberg-size-distributions_',fnames(p,4:7),'-',num2str(berg_cutoff),'mcutoff-L2norm.png'],'png'); 
    %             end
    %         end
            grid on; drawnow;
            % hold off
            
            disp(v1); disp(n1);

        else % if DEM is NOT large enough
            fprintf("==========================================\n")
            fprintf(fnames(p,4:11)); fprintf("\n");
            fprintf("DEM covers too small of an area to get reliable data\n")
            fprintf("==========================================\n")

            % plot just the data
            pl(ref) = loglog(v1,n1,'+','MarkerSize',6,'LineWidth',1,'Color',season_cmap(str2num(fnames(p,8:9)),:)); hold on;
            ylim([0.01,10^8]); %ylim([0.01,10*max(n1)]);
    %         if p == find(yrs == str2num(fnames(p,4:7)),1,'last')
    %             legend(pl,date_names); 
    %             if isinf(norm_type)
    %                 saveas(gcf,[fnames(p,1:2),'-iceberg-size-distributions_',fnames(p,4:7),'-',num2str(berg_cutoff),'mcutoff-INFnorm.png'],'png'); 
    %             else
    %                 saveas(gcf,[fnames(p,1:2),'-iceberg-size-distributions_',fnames(p,4:7),'-',num2str(berg_cutoff),'mcutoff-L2norm.png'],'png'); 
    %             end
    %         end
            grid on; drawnow;

            %fill-in NaNs as outputs
            alpha(p)=NaN;
            c1(p) = NaN; c2(p) = NaN; c3(p) = NaN; c4(p) = NaN;
            data_lims(p,:) = [NaN NaN]; errors(p) = NaN;
        end
        %clear variables
        clearvars -except nfiles berg_cutoff norm_type fitplot time_cmap season_cmap fnames pl ref iceberg_dates date_names a_option yrs alpha c1 c2 c3 c4 data_lims errors;

    end
end
end % close the main function

%% Functions

function n = model(c,v) % Eq. 1 both sides and power law for smallest bergs
% INPUTS:
%   c = vector of parameters
%   v = vector of volume bins
%   k = index to slice for Eq. 1 and the new power law
% OUTPUTS:
%   n = vector of modeled counts in the size distribution

c1 = c(1); c2 = c(2); a = c(3); c3 = c(4); c4 = c(5); c5 = c(6); c6 = c(7); % parse initial guesses
k = c(8); % last "guess" is the k slice index
n = ones(length(v),1)*NaN; % initialize n (same length as v)

if k > 1
    v1 = v(1:k); v2 = v(k+1:end); % slice the volume into two separate vectors
    n(1:k) = c5 * v1.^(-c6); % Power law for smallest bergs
    n(k+1:end) = c1 * v2.^(-a) .* exp(-v2/c2) + c3 * exp(-v2/c4); % Eq. 1
elseif k == 1
    v1 = v(k); v2 = v(k+1:end); % slice v
    n(k) = c5 * v1^(-c6); % Power law for smallest bergs
    n(k+1:end) = c1 * v2.^(-a) .* exp(-v2/c2) + c3 * exp(-v2/c4); % Eq. 1
end
[nanidx] = find(isnan(n)); % find any NaNs
n(nanidx) = 0; % replace NaNs with 0s
end

function [error,c] = fit(v,n,k_cut,i_cut,j_cut,a,norm_type)
% The actual fitting function
% INPUTS:
%   v, n = vectors with the volume bins and counts
%   k_cut = slicing index for fitting power law to smallest bergs
%   i_cut, j_cut = slicing indices for RHS and LHS of Eq. 1
%   a = alpha exponent in Eq. 1 (dictated by fragmentation theory)
%   norm_type = norm to minimize (L2 or maximum)
% OUTPUTS:
%   error = the error associated with the fit(calculated using set norm)
%   c = vector with optimized parameters

addpath('/Users/jukesliu/Documents/BSU_COURSES/2020_Spring/InverseTheory/'); % add inverse theory functions

% Find initial guesses for the optimization:
% 1) new power law fit
vp = v(1:k_cut); np = n(1:k_cut); % slice data

m = length(vp);
A = [ones(m,1),vp]; % initialize matrix A for linearization
pwr = log(np); % linearize np
% b = (A'*A)\(A'*pwr); 
b = pinv(A)*pwr; % Find model using Moore-Penrose Psuedoinverse
cp(1) = exp(b(1)); %c5
cp(2) = b(2); % c6 - power law exponent

% 2) Eq. 1 exponential fit - left side
vl = v(k_cut:i_cut); nl = n(k_cut:i_cut); % slice data
m = length(vl);
A = [ones(m,1),vl];
lhs = log(nl) - log(vl.^(-a)); % linearize by taking the log
% b = (A'*A)\(A'*lhs);
b = pinv(A)*lhs; % psuedoinverse
cl(1) = exp(b(1));
cl(2) = -1/b(2);
cl(3) = a;

% 3) Eq. 1 exponential fit - right side
vr = v(j_cut:end); nr = n(j_cut:end); % slice data vectors
m = length(vr);
A = [ones(m,1),vr];
rhs = log(nr); % linearize by taking the log
% b = (A'*A)\(A'*rhs);
b = pinv(A)*rhs; % psuedoinverse
cr(1) = exp(b(1)); % transform back after initially taking the log
cr(2) = -1/b(2); 

% % now use the linear approach as an initial guess to fit the entire function
c0 = [cl, cr, cp, k_cut]; % Initial guesses
c0 = abs(c0); %making sure initial guess is at least positive in all coeficients

% c = c0; % TOGGLE FOR MANUAL SOLUTIONS VS. LSQCURVEFIT SOLUTIONS
opts = optimset('Display','off');
c = lsqcurvefit(@model, c0, v, n, zeros(length(c0),1), 1e12*ones(length(c0),1),opts); % fitting function with bounds
error = norm((model(c,v)-n)./n,norm_type); % calculate the error associated with the fit

end