%%% Use this code to loop through the iceberg size distribution data that
%%% have fragmentation theory curve parameters (both those with automated
%%% and manual fits) and group data for analysis

%% initialize (run every time)
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/');

%root directory for all glacier folders
root_dir = '/Users/ellynenderlin/Research/NSF_Greenland-Calving/iceberg-calving/';

%automated fit file info
auto_folder = 'models';
auto_filepart = '_parameters';

%manual fit file info
man_folder = 'manually_adjusted_models';
man_filepart = '-parameters-adjusted';

%create a month naming matrix to convert 3-letter months to numeric months
month_array = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};

%cumulative days of year at the start of each month
mo_days = [31 28 31 30 31 30 31 31 30 31 30 31]; cum_days = cumsum(mo_days); cum_days = [0 cum_days(1:11)];
mo_leapdays = [31 29 31 30 31 30 31 31 30 31 30 31]; cum_leapdays = cumsum(mo_leapdays); cum_leapdays = [0 cum_leapdays(1:11)];

%2-letter region flagging (based on alphabetical order of site folders)
region_flag = ['NW'; 'CE'; 'SE'; 'SE'; 'CW'; 'NW'; 'CW'; 'SE'; 'SE'; 'NW'; 'CE'; 'SE'; 'CW'; 'CW'; 'NW'; 'NW'; 'NE'];

%% loop through the files for each site & pull the best fit information

%identify the site folders
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.') && length(sites(i).name) == 2
        sitenames = [sitenames; sites(i).name];
    end
end

%loop through the folders & extract info
disp('Extracting fits...');
for i = 1:length(sitenames)
    cd_to_site = ['cd ',sitenames(i,:),'/']; eval(cd_to_site);
    F(i).site = sitenames(i,:); F(i).region = region_flag(i,:);
    
    %load the average iceberg surface area for each bin
    if i == 1
        distributions = dir('*-iceberg-distribution.txt');
        dist_example = readmatrix(distributions(1).name);
        surfA = dist_example(:,2);
        surfA_binwidth = dist_example(:,3);
    end
    F(i).Abin_mean = surfA'; F(i).Abin_width = surfA_binwidth';
    
    %extract automated fits
    cd_to_folder = ['cd ',auto_folder]; eval(cd_to_folder);
    auto_files = dir(['*',auto_filepart,'*']);
    for j = 1:length(auto_files)
        %identify the month the file was created
        for k = 1:length(month_array)
            if contains(auto_files(j).name(18:20),string(month_array(k,:)))
                file_creation_mo(j) = k;
            end
        end
        
        %identify the day of the year the file was created
        file_creation_doy(j) = str2num(auto_files(j).name(15:16)) + cum_days(file_creation_mo(j));
    end
    [sort_date,sort_ref] = sort(file_creation_doy); %sorted dates of creation
    
    %create a matrix of the fit data, replacing the oldest data with everything more recent
    for j = 1:length(auto_files)
%         T = readtable(auto_files(sort_ref(j)).name); fit_data = table2array(T); 
        fit_data = readmatrix(auto_files(sort_ref(j)).name);
        fit_data(fit_data(:,1)==0,:) = [];
        
        %create a master matrix
        if j == 1; fit_table = fit_data; end
        
        %identify rows of data that need to be replaced with newer estimates
        if j > 1
            for k = 1:size(fit_data,1)
                row_ref = find(fit_table(:,1) - fit_data(k,1) == 0);
                if isempty(row_ref)
                    fit_table = [fit_table; fit_data(k,:)];
                else
                    fit_table(row_ref,:) = fit_data(k,:);
                end
                clear row_ref;
            end
        end
        
        clear T fit_data;
    end
    
    %replace automated fits with manual fits where needed
    cd_to_folder = ['cd ../',man_folder]; eval(cd_to_folder);
    man_files = dir(['*',man_filepart,'.csv']);
    for j = 1:length(man_files)
%         T = readtable(man_files(j).name); fit_data = table2array(T); 
        fit_data = readmatrix(man_files(j).name);
        fit_datestring = man_files(j).name(4:11);
        
        %replace the corresponding automated fit data
        row_ref = find(fit_table(:,1) - str2num(fit_datestring) == 0);
        if size(fit_data,2) == 6
            fit_table(row_ref,2:6) = fit_data(1,2:6); fit_table(row_ref,7:8) = NaN;
        else
            fit_table(row_ref,2:8) = fit_data(1,2:8);
        end
        clear row_ref T fit_data fit_datestring;
    end
    
    %convert all dates to decimal date format
    YYYYMMDD = num2str(fit_table(:,1));
    for j = 1:size(YYYYMMDD,1)
        if mod(str2num(YYYYMMDD(j,1:4)),4) ~= 0
            F(i).dates(j,1) = str2num(YYYYMMDD(j,1:4)) + (cum_days(str2num(YYYYMMDD(j,5:6))) + str2num(YYYYMMDD(j,7:8)))./(max(cum_days)+31);
        else
            F(i).dates(j,1) = str2num(YYYYMMDD(j,1:4)) + (cum_leapdays(str2num(YYYYMMDD(j,5:6))) + str2num(YYYYMMDD(j,7:8)))./(max(cum_days)+31);
        end
    end
    
    %add to over-arching structure
    F(i).alpha = fit_table(:,4); %branching fragmentation exponent
    F(i).c1 = fit_table(:,2); %branching fragmentation constant
    F(i).c2 = fit_table(:,3); %branching fragmentation 'cut-off'
    F(i).c3 = fit_table(:,5); %tabular fragmentation constant
    F(i).c4 = fit_table(:,6); %tabular calvin 'cut-off'
    F(i).c5 = fit_table(:,7); %submarine melt constant
    F(i).beta = fit_table(:,8); %submarine melt exponent
    
    clear file_creation* sort_* fit_table *_files YYYYMMDD;
    cd ../..
    disp(['finished extracting fit time series for ',F(i).site]);
end
save('Greenland-iceberg-fragmentation-curves.mat','F','-v7.3');

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
save('Greenland-iceberg-fragmentation-curves.mat','F','-v7.3');


%% plot curves for each glacier: color-code (1) seasonally & (2) annually
close all;

%load the data
cd /Users/ellynenderlin/Research/NSF_Greenland-Calving/fragmentation-curves
load Greenland-iceberg-fragmentation-curves.mat

%specify years for dataset
start_yr = 2011; end_yr = 2019;

% %specify annual colorramp
% annual_cmap = cmocean('thermal',length([start_yr:1:end_yr]));
% %specify seasonal colorramp
% season_cmap = cmocean('dense',4);

%specify site-specific colorramp
site_cmap = cmocean('curl',length(F));
%specify a regional colorramp
regions = unique(region_flag,'rows');
% region_cmap = cmocean('solar',length(unique(region_flag,'rows')));
region_cmap = [171,217,233; 253,174,97; 44,123,182; 255,255,191; 215,25,28]/255;

%find the intersection of the two fragmentation theory curves & plot
annual_fig = figure; set(annual_fig,'position',[50 50 800 1200]); 
season_fig = figure; set(season_fig,'position',[650 50 800 1200]); 
for i = 1:length(F)
    for j = 1:size(F(i).c1,1)
        branching(j,:) = F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j));
        tabular(j,:) = F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j));
        
        %plot each year as a different color
        figure(annual_fig); 
        subplot(length([start_yr:1:end_yr]),1,floor(F(i).dates(j)-(start_yr-1)));
        if floor(F(i).dates(j)-(start_yr-1)) == length([start_yr:1:end_yr])
           for k = 1:length(regions)
              pl(k) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',region_cmap(k,:),'linewidth',1.5); hold on; 
           end
        end
%         pl(floor(F(i).dates(j)-(start_yr-1))) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',annual_cmap(floor(F(i).dates(j)-(start_yr-1)),:),'linewidth',1.5); hold on;
        loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',region_cmap(strmatch(F(i).region,regions),:),'linewidth',1.5); hold on;
        set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6],'xtick',[100,1000,10000,100000,1000000],'xticklabel',[100,1000,10000,100000,1000000]);
        set(gca,'fontsize',20,'ylim',[10^-13 1]); 
        
        %plot each season as a different color
        figure(season_fig); 
        subplot(4,1,ceil((F(i).dates(j)-floor(F(i).dates(j)))*4));
        if ceil((F(i).dates(j)-floor(F(i).dates(j)))*4) == 4
           for k = 1:length(regions)
              loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color','k','linewidth',2); hold on;
              pr(k) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',region_cmap(k,:),'linewidth',1.5); hold on;
           end
        end
%         pr(ceil((F(i).dates(j)-floor(F(i).dates(j)))*4)) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',season_cmap(ceil((F(i).dates(j)-floor(F(i).dates(j)))*4),:),'linewidth',1.5); hold on;
        loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color','k','linewidth',2); hold on;
        pr(strmatch(F(i).region,regions)) = loglog(F(i).Abin_mean,(branching(j,:)+tabular(j,:))./max(branching(j,:)+tabular(j,:)),'-','color',region_cmap(strmatch(F(i).region,regions),:),'linewidth',1.5); hold on;
        set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6],'xtick',[100,1000,10000,100000,1000000],'xticklabel',[100,1000,10000,100000,1000000]);
        set(gca,'fontsize',20,'ylim',[10^-13 1]); 

        drawnow;
    end
    F(i).branching = branching; F(i).tabular = tabular;
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

%add text labels
figure(annual_fig);
for j = 1:length([start_yr:1:end_yr])
    subplot(length([start_yr:1:end_yr]),1,j);
    text(0.00003*max(get(gca,'xlim')),0.01*max(get(gca,'ylim')),num2str(start_yr-1+j),'fontsize',20);
    if j == length([start_yr:1:end_yr])
    xlabel('Iceberg surface area (m^2)'); ylabel('Normalized count');
plot_pos = get(gca,'position');
leg = legend(pl,regions); set(leg,'location','southoutside','Orientation','horizontal');
set(gca,'position',plot_pos);
    end
end
figure(season_fig); 
subplot(4,1,1);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'a) winter','fontsize',20);
subplot(4,1,2);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'b) spring','fontsize',20);
subplot(4,1,3);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'c) summer','fontsize',20);
subplot(4,1,4);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'d) fall','fontsize',20);
xlabel('Iceberg surface area (km^2)'); ylabel('Normalized count');
plot_pos = get(gca,'position');
legr = legend(pr,regions); set(legr,'location','southoutside','Orientation','horizontal');
set(gca,'position',plot_pos);
clear pl pr;

%save the figures
saveas(annual_fig,'Greenland-annual-iceberg-size-distribution_subplots.png','png');
saveas(season_fig,'Greenland-seasonal-iceberg-size-distribution_subplots.png','png');

%resave with the fragmentation curve data
save('Greenland-iceberg-fragmentation-curves.mat','F','-v7.3');


%% plot quantitative data:(1) fragmentation theory curve intercept, (2) branching fracture cut-off, (3) submarine melt volume loss
close all;

% %load the data
% load Greenland-iceberg-fragmentation-curves.mat

%specify site-specific colorramp
site_cmap = cmocean('curl',length(F));
%specify a regional colorramp
regions = unique(region_flag,'rows');
region_cmap = [171,217,233; 253,174,97; 44,123,182; 255,255,191; 215,25,28]/255;

%find the intersection of the two fragmentation theory curves & plot
figure2 = figure; set(figure2,'position',[50 50 600 300]); 
% sub2a = subplot(2,1,1); sub2b = subplot(2,1,2);
for i = 1:length(F)
    for j = 1:size(F(i).c1,1)
        if ~isempty(polyxpoly(F(i).Abin_mean',F(i).branching(j,:),F(i).Abin_mean',F(i).tabular(j,:)))
            [transition_A,transition_n] = polyxpoly(F(i).Abin_mean',F(i).branching(j,:),F(i).Abin_mean',F(i).tabular(j,:));
            F(i).branchtab_transitionA(j) = transition_A(1); F(i).branchtab_transitionn(j) = transition_n(1); 
        else
            F(i).branchtab_transitionA(j) = NaN; F(i).branchtab_transitionn(j) = NaN;
        end
        clear transition_*;
    end
    
%     %plot each site as a different color
%     subplot(sub2a);
%     pl(i) = plot(F(i).dates,F(i).branchtab_transitionA,'ok','markerfacecolor',site_cmap(i,:),'markersize',12); hold on;
%     %plot each region as a different color
%     subplot(sub2b);
    pr(strmatch(F(i).region,regions)) = plot(F(i).dates,F(i).branchtab_transitionA,'ok','markerfacecolor',region_cmap(strmatch(F(i).region,regions),:),'markersize',12); hold on;
end
% subplot(sub2a);
% leg = legend(pl,sitenames); set(leg,'location','eastoutside');
% yticks = get(gca,'ytick');
% set(gca,'fontsize',20,'yticklabel',yticks./10^6); 
% xlabel('Year'); ylabel('Iceberg surface area (km^2)');
% text(2011.25,0.95*max(get(gca,'ylim')),'a) ','fontsize',20);
% subplot(sub2b);
legr = legend(pr,regions); set(legr,'location','eastoutside');
yticks = get(gca,'ytick');
set(gca,'fontsize',20,'yticklabel',yticks./10^6); 
xlabel('Year'); ylabel('Iceberg surface area (km^2)');
% text(2011.25,0.95*max(get(gca,'ylim')),'b) ','fontsize',20);

%plot the cut-off for the power-law (calving via branching fracture)
figure3 = figure; set(figure3,'position',[450 50 600 300]); clear pl;
% sub3a = subplot(2,1,1); sub3b = subplot(2,1,2);
for i = 1:length(F)
%     %plot each site as a different color
%     subplot(sub3a);
%     pl(i) = plot(F(i).dates,F(i).c2,'ok','markerfacecolor',site_cmap(i,:),'markersize',12); hold on;
%     %plot each region as a different color
%     subplot(sub3b);
    pr(strmatch(F(i).region,regions)) = plot(F(i).dates,F(i).c2,'ok','markerfacecolor',region_cmap(strmatch(F(i).region,regions),:),'markersize',12); hold on;
end
% subplot(sub3a);
% leg = legend(pl,sitenames); set(leg,'location','eastoutside');
% yticks = get(gca,'ytick');
% set(gca,'fontsize',20,'yticklabel',yticks./10^6); 
% xlabel('Year'); ylabel('Branching fracture size cut-off (km^2)');
% text(2011.25,0.95*max(get(gca,'ylim')),'a) ','fontsize',20);
% subplot(sub3b);
legr = legend(pr,regions); set(legr,'location','eastoutside');
yticks = get(gca,'ytick');
set(gca,'fontsize',20,'yticklabel',yticks./10^6); 
xlabel('Year'); ylabel('Branching size cut-off (km^2)');
% text(2011.25,0.95*max(get(gca,'ylim')),'b) ','fontsize',20);
clear pl pr;

%save the figures
saveas(figure2,'Greenland-branching-uncorrelated-calving-transition_scatterplot.png','png');
saveas(figure3,'Greenland-branching-cutoff_scatterplot.png','png');


