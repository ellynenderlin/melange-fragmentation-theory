%%% Use this code to loop through the iceberg size distribution data that
%%% have fragmentation theory curve parameters (both those with automated
%%% and manual fits) and group data for analysis

%% Section 0: Initialize (run every time)
clearvars; close all;
addpath('/Users/icebergs/iceberg-fragmentation/AGU2021');

%specify directories for required files ([root_dir,'/',site_names(i)])
root_dir = '/Users/icebergs/iceberg-fragmentation/'; %include trailing / in file name

%automated fit file info
auto_folder = 'models';
auto_filepart = '_parameters';

%manual fit file info
man_folder = 'manually_adjusted_models';
man_filepart = '-parameters-adjusted';

%create a month naming matrix to convert 3-letter months to numeric months
month_array = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};

%2-letter region flagging (based on alphabetical order of site folders)
site_names = ['AG';'DJ';'HH';'HM';'IB';'IG';'JI';'KB';'KL';'KO';'MD';'MG';'RI';'UM';'UN';'US';'ZI']; %alphabetical site directory list
region_flag = ['NW';'CE';'SE';'NW';'CW';'NW';'CW';'SE';'SE';'NW';'CE';'SE';'CW';'CW';'NW';'NW';'NE'];

%specify a regional colorramp
regions = unique(region_flag,'rows');
region_cmap = [171,217,233; 253,174,97; 44,123,182; 255,255,191; 215,25,28]/255; %same colorblind-friendly color scheme as plot_Greenland_iceberg_calving_fragmentation_fits.m

%% Section 1: Loop through the files for each site & pull the best fit information

 %identify the site folders
%cd(root_dir);
%sites = dir; site_names = [];
 %for i = 1:length(sites)
  %   if ~contains(sites(i).name,'.') && length(sites(i).name) == 2
   %      site_names = [site_names; sites(i).name];
    % end
    %end

%loop through the folders & extract info
disp('Extracting fits...');
for i = 1:length(site_names)
    cd([root_dir,site_names(i,:)]);
    F(i).site = site_names(i,:); F(i).region = region_flag(i,:); %save the site name and region to the data structure
    
    %load the average iceberg surface area for each bin
    if i == 1
        distributions = dir('*-iceberg-distribution.txt');
        dist_example = readmatrix(distributions(1).name); %read the text file as a matrix
        surfA = dist_example(:,2); %grab the surface area (m^2)
        surfA_binwidth = dist_example(:,3); %grab the "bin width" for surface areas
    end
    F(i).Abin_mean = surfA'; F(i).Abin_width = surfA_binwidth'; %save surface are info to the data structure
    
    %extract automated fits
    cd([root_dir,site_names(i,:),'/',auto_folder]);
    auto_files = dir(['*',auto_filepart,'*']); %find the files with the specified name
    %there should only be one file (old ones should now be automatically
    %deleted)
    %create a matrix of the fit data, replacing the oldest data with everything more recent
    %         T = readtable(auto_files(sort_ref(j)).name); fit_data = table2array(T);
    fit_data = readmatrix(auto_files(1).name);
    fit_data(fit_data(:,1)==0,:) = [];
    
    %create a master matrix
%     if j == 1
        fit_table = fit_data; 
%     end
    
    %identify rows of data that need to be replaced with newer estimates
%     if j > 1
%         for k = 1:size(fit_data,1)
%             row_ref = find(fit_table(:,1) - fit_data(k,1) == 0);
%             if isempty(row_ref)
%                 fit_table = [fit_table; fit_data(k,:)];
%             else
%                 fit_table(row_ref,:) = fit_data(k,:);
%             end
%             clear row_ref;
%         end
%     end
%     clear T fit_data;
    
    %replace automated fits with manual fits where needed
    cd([root_dir,site_names(i,:),'/',man_folder]);
    man_files = dir(['*',man_filepart,'.csv']); %find all the manually-adjusted data files
    for j = 1:length(man_files)
%         T = readtable(man_files(j).name); fit_data = table2array(T); 
        fit_data = readmatrix(man_files(j).name); %read in the data as matrices
        fit_datestring = man_files(j).name(4:11); %grab the DEM date from the file name
        
        %replace the corresponding automated fit data
        row_ref = find(fit_table(:,1) - str2num(fit_datestring) == 0); %find the automated fit that needs to be replaced by matching dates
        if ~isempty(row_ref)
            if size(fit_data,2) == 6 %replace fit that does not include submarine melting
                fit_table(row_ref,2:6) = fit_data(1,2:6);
                fit_table(row_ref,7:8) = NaN; %fill in columns 7:8 with NaNs so all dates have the same number of columns
            else %replace fit that includes submarine melting (columns 7:8 have data)
                fit_table(row_ref,2:8) = fit_data(1,2:8);
            end
            %clear row_ref T fit_data fit_datestring;
        else
            error('Date not in the automated fit file. Need to rerun automated fit portion of code then manually check fits.');
        end
    end
    
    %convert all dates to decimal date format
    YYYYMMDD = num2str(fit_table(:,1));
    for j = 1:size(YYYYMMDD,1)
        F(i).dates(j,1) = convert_to_decimaldate(YYYYMMDD(j,:)); %requires convert_to_decimaldate.m function written by Ellyn Enderlin
    end
    
    %add to over-arching structure
    F(i).alpha = fit_table(:,4); %branching fragmentation exponent
    F(i).c1 = fit_table(:,2); %branching fragmentation constant
    F(i).c2 = fit_table(:,3); %branching fragmentation 'cut-off'
    F(i).c3 = fit_table(:,5); %tabular fragmentation constant
    F(i).c4 = fit_table(:,6); %tabular calving 'cut-off'
    F(i).c5 = fit_table(:,7); %submarine melt constant
    F(i).beta = fit_table(:,8); %submarine melt exponent
    
    clear file_creation* sort_* fit_table *_files YYYYMMDD;
    disp(['finished extracting fit time series for ',F(i).site]);
end
save([root_dir,'Greenland-iceberg-fragmentation-curves.mat'],'F','-v7.3');

%filter nonsense results
badcutoff_threshold =  surfA(find(surfA == max(surfA))-4); %semi arbitrary based on typical maximum iceberg sizes around Greenland
for i = 1:length(F)
    if ~isempty(find(F(i).c2>badcutoff_threshold))
%         disp([F(i).site,' (i = ',num2str(i),'), bad ref = ',num2str(find(F(i).c2>badcutoff_threshold))])
        F(i).dates(find(F(i).c2>badcutoff_threshold)) = []; F(i).alpha(find(F(i).c2>badcutoff_threshold)) = [];
        F(i).c1(find(F(i).c2>badcutoff_threshold)) = []; F(i).c3(find(F(i).c2>badcutoff_threshold)) = []; F(i).c4(find(F(i).c2>badcutoff_threshold)) = []; 
        F(i).c5(find(F(i).c2>badcutoff_threshold)) = []; F(i).beta(find(F(i).c2>badcutoff_threshold)) = []; 
        F(i).c2(find(F(i).c2>badcutoff_threshold)) = []; 
    end
end
disp('Filtered bad data based on branching fracture size cut-off');
save([root_dir,'Greenland-iceberg-fragmentation-curves.mat'],'F','-v7.3');


%% Section 2: Plot curves for each glacier: color-code (1) seasonally & (2) annually
close all;

%load the data (if skipping Section 1)
if ~exist('F')
    load([root_dir,'Greenland-iceberg-fragmentation-curves.mat']);
end

%specify years for dataset
start_yr = years(1); end_yr = years(end);

% %specify annual colorramp
% annual_cmap = cmocean('thermal',length([start_yr:1:end_yr]));
% %specify seasonal colorramp
% season_cmap = cmocean('dense',4);

%find the intersection of the two fragmentation theory curves & plot





%Create new variable for 4x4 plots, insert here: (NSEW orientation) Be sure
%to include 2nd direction, just getting it formatted for now: 

%North_fig = figure; (North_fig, 'position', []);
%South_fig = figure; (South_fig, 'position', []);
%East_fig = figure; (East_fig, 'position', []);
%West_fig = figure; (West_fig, 'position', []);
%hold on; 

%for i = 1:length(F) Totally used old code here. Defs feel like an idiot -
%probably not a good plan to reuse code like that.

    %for j = 1:size(F(i).c1,1)
        %calculate the expected iceberg counts for each size bin using the "branching fracture" portion of the fragmentation equation
        %branching(j,:) = F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j)); 
        %calculate the expected iceberg counts for each size bin using the "isolated fracture" portion of the fragmentation equation...
        %produces tabular icebergs
        %tabular(j,:) = F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j));



annual_fig = figure; set(annual_fig,'position',[50 50 800 1200]); 
season_fig = figure; set(season_fig,'position',[650 50 800 1200]); 
for i = 1:length(F)
    for j = 1:size(F(i).c1,1)
        %calculate the expected iceberg counts for each size bin using the "branching fracture" portion of the fragmentation equation
        branching(j,:) = F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j)); 
        %calculate the expected iceberg counts for each size bin using the "isolated fracture" portion of the fragmentation equation...
        %produces tabular icebergs
        tabular(j,:) = F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j));
        
        %create a plot that sorts the data according to year
        figure(annual_fig); 
        subplot(length([start_yr:1:end_yr]),1,floor(F(i).dates(j)-(start_yr-1))); %find which subplot corresponds to the observation year
        %create dummy lines for each region to include in the legend
        if floor(F(i).dates(j)-(start_yr-1)) == length([start_yr:1:end_yr]) 
           for k = 1:length(regions)
              pl(k) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',region_cmap(k,:),'linewidth',1.5); hold on; 
           end
        end
%         pl(floor(F(i).dates(j)-(start_yr-1))) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',annual_cmap(floor(F(i).dates(j)-(start_yr-1)),:),'linewidth',1.5); hold on;
        loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',region_cmap(strmatch(F(i).region,regions),:),'linewidth',1.5); hold on; %plot the fit curves color-coded to region
        set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6],'xtick',[100,1000,10000,100000,1000000],'xticklabel',[100,1000,10000,100000,1000000]); %set the x axis limits and labels
        set(gca,'fontsize',20,'ylim',[10^-13 1]); %set the y axis limits
        
        %create a plot that sorts the data according to season
        figure(season_fig); 
        subplot(4,1,ceil((F(i).dates(j)-floor(F(i).dates(j)))*4)); %find which subplot corresponds to the observation season
        %uncomment if statement below to create dummy lines for each region
        %to include in the legend if the dummy line plotting below the
        %actual data plotting (using loglog) doesn't work
%         if ceil((F(i).dates(j)-floor(F(i).dates(j)))*4) == 4
%            for k = 1:length(regions)
%               loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color','k','linewidth',2); hold on;
%               pr(k) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',region_cmap(k,:),'linewidth',1.5); hold on;
%            end
%         end
%         pr(ceil((F(i).dates(j)-floor(F(i).dates(j)))*4)) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',season_cmap(ceil((F(i).dates(j)-floor(F(i).dates(j)))*4),:),'linewidth',1.5); hold on;
        loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color','k','linewidth',2); hold on; %plot the fit curves color-coded to region
        pr(strmatch(F(i).region,regions)) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),...
            '-','color',region_cmap(strmatch(F(i).region,regions),:),'linewidth',1.5); hold on; %create dummy lines for each region to include in the legend
        set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6],'xtick',[100,1000,10000,100000,1000000],'xticklabel',[100,1000,10000,100000,1000000]); %set the x axis limits and labels
        set(gca,'fontsize',20,'ylim',[10^-13 1]); %set the y axis limits

        drawnow;
    end
    F(i).branching = branching; F(i).tabular = tabular; %add the expected iceberg counts to the data structure
    clear branching tabular;
end
% subplot(sub1a);
% leg = legend(pl,{num2str([2011:1:2019]')}); set(leg,'location','eastoutside');
% set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6]); xticks = get(gca,'xtick');
% set(gca,'fontsize',20,'xticklabel',xticks./10^6,'ylim',[10^-13 1]); 
% xlabel('Iceberg surface area (km^2)'); ylabel('Normalized count');
% text(0.00003*max(get(gca,'xlim')),0.4*max(get(gca,'ylim')),'a) ','fontsize',20);
% subplot(sub1b);
% legr = legend(pr,'winter','spring','summer','fall'); set(legr,'location','eastoutside');
% set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6]); xticks = get(gca,'xtick');
% set(gca,'fontsize',20,'xticklabel',xticks./10^6,'ylim',[10^-13 1]); 
% xlabel('Iceberg surface area (km^2)'); ylabel('Normalized count');
% text(0.00003*max(get(gca,'xlim')),0.4*max(get(gca,'ylim')),'b) ','fontsize',20);

%format the annual plot figure
figure(annual_fig);
for j = 1:length([start_yr:1:end_yr])
    subplot(length([start_yr:1:end_yr]),1,j);
    text(0.00003*max(get(gca,'xlim')),0.01*max(get(gca,'ylim')),num2str(start_yr-1+j),'fontsize',20);
    %only add x and y axis labels to the lower left plot (may need to
    %change if statement to length([start_yr:1:end_yr])-1)
    if j == length([start_yr:1:end_yr])
        xlabel('Iceberg surface area (m^2)'); ylabel('Normalized count');
        plot_pos = get(gca,'position');
        leg = legend(pl,regions); set(leg,'location','southoutside','Orientation','horizontal'); %move and reorient the legend
        set(gca,'position',plot_pos); %fix the plot position so it is back in the default pre-legend location
    end
end

%add season labels to the seasonal subplots
figure(season_fig); 
subplot(4,1,1);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'a) winter','fontsize',20);
subplot(4,1,2);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'b) spring','fontsize',20);
subplot(4,1,3);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'c) summer','fontsize',20);
subplot(4,1,4);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'d) fall','fontsize',20);
%format the seasonal plot figure
subplot(4,1,3); xlabel('Iceberg surface area (km^2)'); ylabel('Normalized count'); %add axis labels to the lower left panel
plot_pos = get(gca,'position');
legr = legend(pr,regions); %add a legend to distinguish region colors
set(legr,'location','southoutside','Orientation','horizontal'); %move and reorient the legend
set(gca,'position',plot_pos); %fix the plot position so it is back in the default pre-legend location
clear pl pr;

%save the figures
saveas(annual_fig,[root_dir,'Greenland-annual-iceberg-size-distribution_subplots.eps'],'epsc');
saveas(annual_fig,[root_dir,'Greenland-annual-iceberg-size-distribution_subplots.png'],'png');
saveas(season_fig,[root_dir,'Greenland-seasonal-iceberg-size-distribution_subplots.eps'],'epsc');
saveas(season_fig,[root_dir,'Greenland-seasonal-iceberg-size-distribution_subplots.png'],'png');

%resave with the fragmentation curve data
save([root_dir,'Greenland-iceberg-fragmentation-curves.mat'],'F','-v7.3');


%% Section 3: Plot quantitative data:(1) fragmentation theory curve intercept, (2) branching fracture cut-off, (3) submarine melt volume loss
close all;

%load the data (if skipping Sections 1-2)
if ~exist('F')
    load([root_dir,'Greenland-iceberg-fragmentation-curves.mat']);
end

%specify site-specific colorramp
site_cmap = cmocean('curl',length(F));

%find the intersection of the two fragmentation theory curves & plot: this
%tells us the size threshold where iceberg sizes switch from being
%dominated by branching fracture processes to isolated fracture
figure2 = figure; set(figure2,'position',[50 50 600 300]); 
for i = 1:length(F)
    for j = 1:size(F(i).c1,1)
        if ~isempty(polyxpoly(F(i).Abin_mean',F(i).branching(j,:),F(i).Abin_mean',F(i).tabular(j,:))) %if the curves intersect, solve for the intersection
            [transition_A,transition_n] = polyxpoly(F(i).Abin_mean',F(i).branching(j,:),F(i).Abin_mean',F(i).tabular(j,:)); %use polyxpoly build-in function to solve for curve intersection
            F(i).branchtab_transitionA(j) = transition_A(1); F(i).branchtab_transitionn(j) = transition_n(1); %add the area and normalized iceberg count for the curve intersection to the data structure
        else
            F(i).branchtab_transitionA(j) = NaN; F(i).branchtab_transitionn(j) = NaN; %add NaNs for the area and normalized iceberg count for the curve intersection to the data structure
        end
        clear transition_*;
    end
    %plot a timeseries of iceberg area size for the intersection of the curves, with marker faces colored by region
    pr(strmatch(F(i).region,regions)) = plot(F(i).dates,F(i).branchtab_transitionA,'ok','markerfacecolor',region_cmap(strmatch(F(i).region,regions),:),'markersize',12); hold on;
end
%format the plot
legr = legend(pr,regions); set(legr,'location','eastoutside');
yticks = get(gca,'ytick');
set(gca,'fontsize',20,'yticklabel',yticks./10^6); 
xlabel('Year'); ylabel('Iceberg surface area (km^2)');


%plot the cut-off for the power-law (calving via branching fracture): this
%tells us the threshold for iceberg area above which the iceberg size
%distribution no longer adheres to a power law distribution
figure3 = figure; set(figure3,'position',[450 50 600 300]); clear pl;
for i = 1:length(F)
    pr(strmatch(F(i).region,regions)) = plot(F(i).dates,F(i).c2,'ok','markerfacecolor',region_cmap(strmatch(F(i).region,regions),:),'markersize',12); hold on;
end
legr = legend(pr,regions); set(legr,'location','eastoutside');
yticks = get(gca,'ytick');
set(gca,'fontsize',20,'yticklabel',yticks./10^6); 
xlabel('Year'); ylabel('Branching size cut-off (km^2)');
clear pl pr;

%save the figures
saveas(figure2,[root_dir,'Greenland-branching-uncorrelated-calving-transition_scatterplot.eps'],'epsc');
saveas(figure2,[root_dir,'Greenland-branching-uncorrelated-calving-transition_scatterplot.png'],'png');
saveas(figure3,[root_dir,'Greenland-branching-cutoff_scatterplot.eps'],'epsc');
saveas(figure3,[root_dir,'Greenland-branching-cutoff_scatterplot.png'],'png');


