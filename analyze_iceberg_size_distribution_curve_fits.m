%%% Use this code to loop through the iceberg size distribution data that
%%% have fragmentation theory curve parameters (both those with automated
%%% and manual fits) and group data for analysis

%Section 0: Initializes the code. This should be the only subsection that
%needs to be modified by the user. It is set-up to pull data from the
%automatically-created "models" and "manually_adjusted_models" directories
%within each site's directory structure, with 16 sites used for code
%generation and analysis by the developer (see site_abbrevs variable
%specifying site directory names and site_names variable for their official
%names). Relies on independent cmocean and ArcticMappingTools that must be
%added to the Matlab code path.

%Section 1: Creates a matfile containing all observed size distributions
%and best model fits (manual overwrite automated) for all sites saved in a single structure, F, called
%Greenland-iceberg-fragmentation-curves.mat.

%Sections 2-6: Visualize the data in different ways, pulling data from the
%compiled site structure, F, created in Section 1.

%% Section 0: Initialize (run every time)
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/');

%specify directories for required files ([root_dir,'/',site_abbrevs(i)])
% root_dir = '/Volumes/CALVING/Greenland_icebergs/iceberg-fragmentation/'; %include trailing / in file name
root_dir = '/Users/ellynenderlin/Research/NSF_Greenland-Calving/iceberg-calving/';
cd(root_dir);

%automated fit file info
auto_folder = 'models';
auto_filepart = '_parameters';

%manual fit file info
man_folder = 'manually_adjusted_models';
man_filepart = '-parameters-adjusted';

%create a month naming matrix to convert 3-letter months to numeric months
month_array = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};
years = 2011:1:2020; start_yr = years(1); end_yr = years(end);

%2-letter region flagging (based on alphabetical order of site folders)
site_abbrevs = ['AG';'HH';'HM';'IB';'IG';'JI';'KB';'KL';'KO';'MD';'MG';'RI';'UM';'UN';'US';'ZI']; %alphabetical site directory list
site_names = [{'Alison'},{'Helheim'},{'Ullip Sermia'},{'Salliarutsip Sermia'},{'Illulip Sermia'},...
    {'Sermeq Kujalleq'},{'Koge Bugt'},{'Kangerlussuaq'},{'Kong Oscar'},{'Magga Dan'},...
    {'Nigertiip Apusiia'},{'Kangilliup Sermia'},{'Umiammakku Sermia'},{'Upernavik North'},{'Upernavik South'},{'Zachariae Isstrom'}];
region_flag = ['NW';'SE';'NW';'CW';'NW';'CW';'SE';'SE';'NW';'CE';'SE';'CW';'CW';'NW';'NW';'NE'];
site_geog_order = [3;12;1;7;4;10;11;14;2;15;13;9;8;5;6;16]; %geographic order, counterclockwise from NW
%NOTE: data are added to F structure according to index specified by
%site_geog_order and plot locations and letters below are for those sorted data
%create a vector specifying subplot locations according to geography
site_plotlocs = ones(size(site_geog_order));
site_plotlocs(1:10) = 2*[1:1:10]'-1; %based on 10 sites along W coast
site_plotlocs(11:end-1) = [11:1:length(site_plotlocs)-1]' - 3*([11:1:length(site_plotlocs)-1]'-14);
site_plotlocs(end) = 2;
%letters for labels
site_plotletts = [{'a)'},{'b)'},{'c)'},{'d)'},{'e)'},{'f)'},{'g)'},{'h)'},{'i)'},{'j)'},{'k)'},{'l)'},{'m)'},{'n)'},{'o)'},{'p)'}];

%specify a regional colorramp
regions = unique(region_flag,'rows');
region_cmap = [171,217,233; 253,174,97; 44,123,182; 255,255,191; 215,25,28]/255; %same colorblind-friendly color scheme as plot_Greenland_iceberg_calving_fragmentation_fits.m

%specify annual colorramp
annual_cmap = cmocean('thermal',length([start_yr:1:end_yr]));
%specify seasonal colorramp
season_cmap = cmocean('curl',6); season_cmap = season_cmap(2:end-1,:);
%specify daily colorramp
% day_cmap = cmocean('phase',365);
accum_cmap = cmocean('ice',233); ablat_cmap = flipud(cmocean('solar',162));
day_cmap = [accum_cmap(103:end-10,:); ablat_cmap(11:end,:); accum_cmap(11:102,:)];


disp('Initialized iceberg fragmentation data compilation code');
%% Section 1: Loop through the files for each site & pull the best fit information
disp('Compile observed & fit iceberg size distributions');

 %identify the site folders
%cd(root_dir);
%sites = dir; site_abbrevs = [];
 %for i = 1:length(sites)
  %   if ~contains(sites(i).name,'.') && length(sites(i).name) == 2
   %      site_abbrevs = [site_abbrevs; sites(i).name];
    % end
    %end

%loop through the folders & extract info
disp('Extracting observations & fits...');
for i = 1:max(site_geog_order)
    cd([root_dir,site_abbrevs(i,:)]);
    F(site_geog_order(i)).site_name = site_names(:,i); 
    F(site_geog_order(i)).site_abbrev = site_abbrevs(i,:); F(site_geog_order(i)).region = region_flag(i,:); %save the site name and region to the data structure
    
    %load the center coordinates for the conservative melange polygon
    load([site_abbrevs(i,:),'-melange-masks.mat']);
    F(site_geog_order(i)).x = nanmean(melmask.uncropped.x); F(site_geog_order(i)).y = nanmean(melmask.uncropped.y); 
    clear melmask;
    
    %load the average iceberg surface area for each bin
    distributions = dir('*-iceberg-distribution.txt');
    for j = 1:length(distributions)
        dist_example = readmatrix(distributions(j).name); %read the text file as a matrix
        bin_no(:,j) = dist_example(:,1); %grab the number of icebergs for the surface area bin
        surfA = dist_example(:,2); %grab the surface area (m^2)
        surfA_binwidth = dist_example(:,3); %grab the "bin width" for surface areas
    end
    F(site_geog_order(i)).Abin_mean = surfA'; F(site_geog_order(i)).Abin_width = surfA_binwidth'; %save surface are info to the data structure
    
    %compile observed numbers of icebergs for each bin
    F(site_geog_order(i)).Abin_no = bin_no';
    
    %extract automated fits
    cd([root_dir,site_abbrevs(i,:),'/',auto_folder]);
    auto_files = dir(['*',auto_filepart,'*']); %find the files with the specified name
    fit_data = readmatrix(auto_files(1).name);
    fit_data(fit_data(:,1)==0,:) = [];
    
    %create a master matrix
    fit_table = fit_data; 
    
    %replace automated fits with manual fits where needed
    cd([root_dir,site_abbrevs(i,:),'/',man_folder]);
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
        F(site_geog_order(i)).dates(j,1) = convert_to_decimaldate(YYYYMMDD(j,:)); %requires convert_to_decimaldate.m function written by Ellyn Enderlin
    end
    
    %add to over-arching structure
    F(site_geog_order(i)).alpha = fit_table(:,4); %branching fragmentation exponent
    F(site_geog_order(i)).c1 = fit_table(:,2); %branching fragmentation constant
    F(site_geog_order(i)).c2 = fit_table(:,3); %branching fragmentation 'cut-off'
    F(site_geog_order(i)).c3 = fit_table(:,5); %tabular fragmentation constant
    F(site_geog_order(i)).c4 = fit_table(:,6); %tabular calving 'cut-off'
    F(site_geog_order(i)).c5 = fit_table(:,7); %submarine melt constant
    F(site_geog_order(i)).beta = fit_table(:,8); %submarine melt exponent
    
    clear bin_no file_creation* sort_* fit_table *_files YYYYMMDD;
    disp(['finished extracting fit time series for ',F(site_geog_order(i)).site_abbrev]);
end
save([root_dir,'Greenland-iceberg-fragmentation-curves.mat'],'F','-v7.3');

%filter nonsense results
badcutoff_threshold =  surfA(find(surfA == max(surfA))-4); %semi arbitrary based on typical maximum iceberg sizes around Greenland
for i = 1:length(F)
    if ~isempty(find(F(site_geog_order(i)).c2>badcutoff_threshold))
%         disp([F(site_geog_order(i)).site_abbrev,' (i = ',num2str(i),'), bad ref = ',num2str(find(F(site_geog_order(i)).c2>badcutoff_threshold))])
        F(site_geog_order(i)).dates(find(F(site_geog_order(i)).c2>badcutoff_threshold)) = []; F(site_geog_order(i)).alpha(find(F(site_geog_order(i)).c2>badcutoff_threshold)) = [];
        F(site_geog_order(i)).c1(find(F(site_geog_order(i)).c2>badcutoff_threshold)) = []; F(site_geog_order(i)).c3(find(F(site_geog_order(i)).c2>badcutoff_threshold)) = []; F(site_geog_order(i)).c4(find(F(site_geog_order(i)).c2>badcutoff_threshold)) = []; 
        F(site_geog_order(i)).c5(find(F(site_geog_order(i)).c2>badcutoff_threshold)) = []; F(site_geog_order(i)).beta(find(F(site_geog_order(i)).c2>badcutoff_threshold)) = []; 
        F(site_geog_order(i)).c2(find(F(site_geog_order(i)).c2>badcutoff_threshold)) = []; 
    end
end
disp('Filtered bad data based on branching fracture size cut-off');
save([root_dir,'Greenland-iceberg-fragmentation-curves.mat'],'F','-v7.3');
disp('Saved compiled size distribution data to a single mat file');

%% Section 2: Plot observed size distributions for each glacier
close all;

%load the data (if skipping Section 1)
if ~exist('F')
    load([root_dir,'Greenland-iceberg-fragmentation-curves.mat']);
end


%create figures showing observations for each site
for i = 1:length(F)
    %create the template figures
    annual_sitefig = figure; set(annual_sitefig,'position',[50 50 800 600]);
    for k = 1:size(annual_cmap,1)
        dpa(k) = loglog(F(i).Abin_mean(1:2),F(i).Abin_no(1,1:2),'-','linewidth',1.5,'color',annual_cmap(k,:)); hold on;
    end
    season_sitefig = figure; set(season_sitefig,'position',[50 650 800 600]);
    for k = 1:size(season_cmap,1)
        dps(k) = loglog(F(i).Abin_mean(1:2),F(i).Abin_no(1,1:2),'-','linewidth',1.5,'color',season_cmap(k,:)); hold on;
    end
    
    %plot the observed size distributions
    for j = 1:size(F(i).Abin_no,1)
        figure(annual_sitefig);
        loglog(F(i).Abin_mean,F(i).Abin_no(j,:),'-','linewidth',1.5,'color',annual_cmap(floor(F(i).dates(j)-start_yr)+1,:)); hold on;
        figure(season_sitefig);
        loglog(F(i).Abin_mean,F(i).Abin_no(j,:),'-','linewidth',1.5,'color',season_cmap(ceil(4*(F(i).dates(j)-floor(F(i).dates(j)))),:)); hold on;
    end
    %add labels to the annual profiles
    figure(annual_sitefig); grid on; ylims = get(gca,'ylim');
    set(gca,'xlim',[0 2e6],'xtick',[100,1000,10000,100000,1000000]); %set the x axis limits and labels
    set(gca,'fontsize',16,'ylim',[0.1 max(ylims)]); %set the y axis limits
    xlabel('Iceberg surface area (m^2)'); ylabel('Count');
    lega = legend(dpa,num2str([start_yr:1:end_yr]'));
    drawnow;
    %add labels to the seasonal profiles
    figure(season_sitefig); grid on; ylims = get(gca,'ylim');
    set(gca,'xlim',[0 2e6],'xtick',[100,1000,10000,100000,1000000]); %set the x axis limits and labels
    set(gca,'fontsize',16,'ylim',[0.1 max(ylims)]); %set the y axis limits
    xlabel('Iceberg surface area (m^2)'); ylabel('Count');
    legs = legend(dps,'winter','spring','summer','autumn');
    drawnow;
    
    %save the figures then close
    saveas(annual_sitefig,[root_dir,F(i).site_abbrev,'/',F(i).site_abbrev,'-annual-iceberg-size-distribution_plot.png'],'png');
    saveas(season_sitefig,[root_dir,F(i).site_abbrev,'/',F(i).site_abbrev,'-seasonal-iceberg-size-distribution_plot.png'],'png');
    close(annual_sitefig); close(season_sitefig);
    clear dpa dps lega legs;
end

%create figures showing fitted models for each site
for i = 1:length(F)
    %create the template figures
    annual_sitefig = figure; set(annual_sitefig,'position',[50 50 800 600]);
    for k = 1:size(annual_cmap,1)
        dpa(k) = loglog(F(i).Abin_mean(1:2),F(i).c1(1).*F(i).Abin_mean(1:2).^-F(i).alpha(1).*exp(-F(i).Abin_mean(1:2)./F(i).c2(1))+F(i).c3(1).*exp(-F(i).Abin_mean(1:2)./F(i).c4(1)),'-','linewidth',1.5,'color',annual_cmap(k,:)); hold on;
    end
    season_sitefig = figure; set(season_sitefig,'position',[50 650 800 600]);
    for k = 1:size(season_cmap,1)
        dps(k) = loglog(F(i).Abin_mean(1:2),F(i).c1(1).*F(i).Abin_mean(1:2).^-F(i).alpha(1).*exp(-F(i).Abin_mean(1:2)./F(i).c2(1))+F(i).c3(1).*exp(-F(i).Abin_mean(1:2)./F(i).c4(1)),'-','linewidth',1.5,'color',season_cmap(k,:)); hold on;
    end
    
    %plot the observed size distributions
    for j = 1:size(F(i).Abin_no,1)
        figure(annual_sitefig);
        loglog(F(i).Abin_mean,F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j))+F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j)),'-','linewidth',1.5,'color',annual_cmap(floor(F(i).dates(j)-start_yr)+1,:)); hold on;
        figure(season_sitefig);
        loglog(F(i).Abin_mean,F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j))+F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j)),'-','linewidth',1.5,'color',season_cmap(ceil(4*(F(i).dates(j)-floor(F(i).dates(j)))),:)); hold on;
    end
    %add labels to the annual profiles
    figure(annual_sitefig); grid on; ylims = get(gca,'ylim');
    set(gca,'xlim',[0 2e6],'xtick',[100,1000,10000,100000,1000000]); %set the x axis limits and labels
    set(gca,'fontsize',16,'ylim',[10^-6 10^3]); %set the y axis limits
    xlabel('Iceberg surface area (m^2)'); ylabel('Count');
    lega = legend(dpa,num2str([start_yr:1:end_yr]'));
    drawnow;
    %add labels to the seasonal profiles
    figure(season_sitefig); grid on; ylims = get(gca,'ylim');
    set(gca,'xlim',[0 2e6],'xtick',[100,1000,10000,100000,1000000]); %set the x axis limits and labels
    set(gca,'fontsize',16,'ylim',[10^-6 10^3]); %set the y axis limits
    xlabel('Iceberg surface area (m^2)'); ylabel('Count');
    legs = legend(dps,'winter','spring','summer','autumn');
    drawnow;
    
    %save the figures then close
    saveas(annual_sitefig,[root_dir,F(i).site_abbrev,'/',F(i).site_abbrev,'-annual-iceberg-size-fragmentation_plot.png'],'png');
    saveas(season_sitefig,[root_dir,F(i).site_abbrev,'/',F(i).site_abbrev,'-seasonal-iceberg-size-fragmentation_plot.png'],'png');
    close(annual_sitefig); close(season_sitefig);
    clear dpa dps lega legs;
end
disp('Created annual & seasonal plots of observations and fragmentation fits for each site');


%% Section 3: Plot lumped data normalized several ways
close all;

%sort by accumulation or ablation season but without any normalization
season_fig = figure; set(season_fig,'position',[50 50 800 400]);
subacc = subplot(1,2,1); subabl = subplot(1,2,2);
%add the data to the figure template
for i = 1:length(F)
    %plot the observed size distributions
    for j = 1:size(F(i).Abin_no,1)
        if round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) <= 121 || round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) > 274
            subplot(subacc);
        else
            subplot(subabl);
        end
        loglog(F(i).Abin_mean,F(i).Abin_no(j,:),'-','linewidth',1.5,'color',day_cmap(round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))),:)); hold on;
        drawnow;
    end

end
%format subplots
subplot(subacc); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.1*pos(3) 1.05*pos(4)]);
grid on; ylims = get(gca,'ylim');
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[0.1 max(ylims)],'ytick',[1e0 1e3 1e6]); %set the y axis limits
text(70,0.2,'a) accumulation season','FontSize',16);
xlabel('Iceberg surface area (m^2)'); ylabel('Iceberg count');
subplot(subabl); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.1*pos(3) 1.05*pos(4)]);
grid on; ylims = get(gca,'ylim');
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[0.1 max(ylims)],'ytick',[1e0 1e3 1e6]); %set the y axis limits
set(gca,'yticklabel',[]);
text(70,0.2,'b) ablation season','FontSize',16);
xlabel('Iceberg surface area (m^2)'); 
drawnow;
clear sub*;
saveas(season_fig,[root_dir,'Greenland-daily-iceberg-size-fragmentation_lumped-plot.png'],'png');

%normalize by total area (instead of count, the number is iceberg concentration (count/km^2))
seasoncon_fig = figure; set(seasoncon_fig,'position',[550 50 800 900]);
subaccW = subplot(2,2,1); subablW = subplot(2,2,3);
subaccE = subplot(2,2,2); subablE = subplot(2,2,4);
%add the data to the figure template
for i = 1:length(F)
    %plot the observed size distributions
    for j = 1:size(F(i).Abin_no,1)
        if round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) <= 121 || round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) > 274
            if i <=10; subplot(subaccW); else subplot(subaccE); end
        else
            if i <=10; subplot(subablW); else subplot(subablE); end
        end
        loglog(F(i).Abin_mean,F(i).Abin_no(j,:)./(sum(F(i).Abin_mean.*F(i).Abin_no(j,:))/10^6),'-','linewidth',1.5,'color',day_cmap(round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))),:)); hold on;
        drawnow;
    end

end
%format subplots
subplot(subablW); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
grid on; 
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[5e-4 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.8e-3,'b) west: ablation season','FontSize',16);
ylabel('Iceberg concentration (count/km^2)'); xlabel('Iceberg surface area (m^2)'); 
subplot(subablE); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
grid on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[5e-4 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
set(gca,'yticklabel',[]);
text(70,0.8e-3,'d) east: ablation season','FontSize',16);
xlabel('Iceberg surface area (m^2)'); 
subplot(subaccW); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
grid on; 
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[5e-4 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.8e-3,'a) west: accumulation season','FontSize',16);
ylabel('Iceberg concentration (count/km^2)'); 
subplot(subaccE); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
grid on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[5e-4 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.8e-3,'c) east: accumulation season','FontSize',16);
set(gca,'yticklabel',[]);
drawnow;
clear sub*;
saveas(seasoncon_fig,[root_dir,'Greenland-daily-iceberg-size-concentration_lumped-plot.png'],'png');


%normalize count by bin width
seasonnorm_fig = figure; set(seasonnorm_fig,'position',[1050 50 800 900]);
subaccW = subplot(2,2,1); subablW = subplot(2,2,3);
subaccE = subplot(2,2,2); subablE = subplot(2,2,4);
%add the data to the figure template
for i = 1:length(F)
    %plot the observed size distributions
    for j = 1:size(F(i).Abin_no,1)
        if round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) <= 121 || round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) > 274
            if i <=10; subplot(subaccW); else subplot(subaccE); end
        else
            if i <=10; subplot(subablW); else subplot(subablE); end
        end
        loglog(F(i).Abin_mean,F(i).Abin_no(j,:)./F(i).Abin_width,'-','linewidth',1.5,'color',day_cmap(round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))),:)); hold on;
        drawnow;
    end

end
%format subplots
subplot(subablW); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
grid on; 
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'b) west: ablation season','FontSize',16);
ylabel('Iceberg number density'); xlabel('Iceberg surface area (m^2)'); 
subplot(subablE); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
grid on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
set(gca,'yticklabel',[]);
text(70,0.2e-4,'d) east: ablation season','FontSize',16);
xlabel('Iceberg surface area (m^2)'); 
subplot(subaccW); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
grid on; 
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'a) west: accumulation season','FontSize',16);
ylabel('Iceberg number density'); 
subplot(subaccE); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
grid on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'c) east: accumulation season','FontSize',16);
set(gca,'yticklabel',[]);
drawnow;
clear sub*;
saveas(seasonnorm_fig,[root_dir,'Greenland-daily-iceberg-size-density_lumped-plot.png'],'png');


%plot each glacier's mean normalized profile for each season
meancon_fig = figure; set(meancon_fig,'position',[550 550 800 900]);
subaccWc = subplot(2,2,1); subablWc = subplot(2,2,3); subaccEc = subplot(2,2,2); subablEc = subplot(2,2,4);
meannorm_fig = figure; set(meannorm_fig,'position',[1050 550 800 900]);
subaccWn = subplot(2,2,1); subablWn = subplot(2,2,3); subaccEn = subplot(2,2,2); subablEn = subplot(2,2,4);
%set up dummy matrices for mean profiles
accumW_no = []; accumW_con = []; accumE_no = []; accumE_con = [];
ablatW_no = []; ablatW_con = []; ablatE_no = []; ablatE_con = [];
%add the data to the figure template
for i = 1:length(F)
    accum_Abin = F(i).Abin_mean; ablat_Abin = F(i).Abin_mean; 
    accum_no = []; ablat_no = [];
    accum_con = []; ablat_con = [];
    
    %parse the data by season
    for j = 1:size(F(i).Abin_no,1)
        if round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) <= 121 || round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) > 274
            accum_no = [accum_no; F(i).Abin_no(j,:)];
            accum_con = [accum_con; F(i).Abin_no(j,:)./(sum(accum_Abin.*F(i).Abin_no(j,:))/10^6)];
        else
            ablat_no = [ablat_no; F(i).Abin_no(j,:)];
            ablat_con = [ablat_con; F(i).Abin_no(j,:)./(sum(ablat_Abin.*F(i).Abin_no(j,:))/10^6)];
        end

    end
    
    %plot the mean profiles
    figure(meancon_fig); %iceberg concentration profiles
    if i <=10 %accumulation season
        subplot(subaccWc); accumW_con = [accumW_con; nanmean(accum_con,1)]; 
    else
        subplot(subaccEc); accumE_con = [accumE_con; nanmean(accum_con,1)]; 
    end 
    loglog(accum_Abin,accum_con,'-','linewidth',1,'color',[0.75 0.75 0.75]); hold on;
    if i <=10 %ablation season
        subplot(subablWc); ablatW_con = [ablatW_con; nanmean(ablat_con,1)]; 
    else
        subplot(subablEc); ablatE_con = [ablatE_con; nanmean(ablat_con,1)]; 
    end 
    loglog(ablat_Abin,ablat_con,'-','linewidth',1,'color',[0.75 0.75 0.75]); hold on;
    figure(meannorm_fig); %iceberg 'density'... count normalized by bin size
    if i <=10 %accumulation season
        subplot(subaccWn); accumW_no = [accumW_no; nanmean(accum_no,1)]; 
    else
        subplot(subaccEn); accumE_no = [accumE_no; nanmean(accum_no,1)]; 
    end 
    loglog(accum_Abin,accum_no./F(i).Abin_width,'-','linewidth',1,'color',[0.75 0.75 0.75]); hold on;
    if i <=10 %ablation season
        subplot(subablWn); ablatW_no = [ablatW_no; nanmean(ablat_no,1)]; 
    else
        subplot(subablEn); ablatE_no = [ablatE_no; nanmean(ablat_no,1)]; 
    end 
    loglog(ablat_Abin,ablat_no./F(i).Abin_width,'-','linewidth',1,'color',[0.75 0.75 0.75]); hold on;
    drawnow;
    
    %clear variables
    clear accum_Abin ablat_Abin accum_no accum_con ablat_no ablat_con;
end
%format subplots
figure(meancon_fig); 
subplot(subablWc); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,ablatW_con,'-','linewidth',2,'color',nanmean(ablat_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(ablatW_con,1),'-k','linewidth',2); hold on; 
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[5e-4 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.8e-3,'b) west: ablation season','FontSize',16);
ylabel('Iceberg concentration (count/km^2)'); xlabel('Iceberg surface area (m^2)'); 
subplot(subablEc); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,ablatE_con,'-','linewidth',2,'color',nanmean(ablat_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(ablatE_con,1),'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[5e-4 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.8e-3,'d) east: ablation season','FontSize',16);
set(gca,'yticklabel',[]);
xlabel('Iceberg surface area (m^2)'); 
subplot(subaccWc); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,accumW_con,'-','linewidth',2,'color',nanmean(accum_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(accumW_con,1),'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[5e-4 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.8e-3,'a) west: accumulation season','FontSize',16);
ylabel('Iceberg concentration (count/km^2)'); 
subplot(subaccEc); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,accumE_con,'-','linewidth',2,'color',nanmean(accum_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(accumE_con,1),'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[5e-4 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.8e-3,'c) east: accumulation season','FontSize',16);
set(gca,'yticklabel',[]);
drawnow;
figure(meannorm_fig); 
subplot(subablWn); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,ablatW_no./F(i).Abin_width,'-','linewidth',2,'color',nanmean(ablat_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(ablatW_no,1)./F(i).Abin_width,'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'b) west: ablation season','FontSize',16);
ylabel('Iceberg number density'); xlabel('Iceberg surface area (m^2)'); 
subplot(subablEn); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,ablatE_no./F(i).Abin_width,'-','linewidth',2,'color',nanmean(ablat_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(ablatE_no,1)./F(i).Abin_width,'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'d) east: ablation season','FontSize',16);
set(gca,'yticklabel',[]);
xlabel('Iceberg surface area (m^2)'); 
subplot(subaccWn); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,accumW_no./F(i).Abin_width,'-','linewidth',2,'color',nanmean(accum_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(accumW_no,1)./F(i).Abin_width,'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'a) west: accumulation season','FontSize',16);
ylabel('Iceberg number density'); 
subplot(subaccEn); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,accumE_no./F(i).Abin_width,'-','linewidth',2,'color',nanmean(accum_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(accumE_no,1)./F(i).Abin_width,'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'c) east: accumulation season','FontSize',16);
set(gca,'yticklabel',[]);
drawnow;
saveas(meancon_fig,[root_dir,'Greenland-mean-seasonal-iceberg-size-concentration_lumped-plot.png'],'png');
saveas(meancon_fig,[root_dir,'Greenland-mean-seasonal-iceberg-size-concentration_lumped-plot.eps'],'epsc');
saveas(meannorm_fig,[root_dir,'Greenland-mean-seasonal-iceberg-size-density_lumped-plot.png'],'png');

%recreate the same normalized lumped plot but for model fits
modelnorm_fig = figure; set(modelnorm_fig,'position',[1550 550 800 900]);
subaccWd = subplot(2,2,1); subablWd = subplot(2,2,3); subaccEd = subplot(2,2,2); subablEd = subplot(2,2,4);
%set up dummy matrices for mean profiles
accumW_dens = []; accumE_dens = [];
ablatW_dens = []; ablatE_dens = [];
%add the data to the figure template
for i = 1:length(F)
    accum_Abin = F(i).Abin_mean; ablat_Abin = F(i).Abin_mean; 
    accum_dens = []; ablat_dens = [];
    
    %parse the data by season
    for j = 1:size(F(i).Abin_no,1)
        if round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) <= 121 || round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) > 274
            accum_dens = [accum_dens; (F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j)))+(F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j)))];
        else
            ablat_dens = [ablat_dens; (F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j)))+(F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j)))];
        end
        
        %flag really bad profiles
        if max((F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j)))+(F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j)))) < 10^1;
            disp(['Bad modeled fit for ',F(i).site_abbrev,' for date ',num2str(F(i).dates(j))]);
        end
    end
    
    
    %plot the mean profiles of modeled iceberg 'density'... count normalized by bin size
    if i <=10 %accumulation season
        subplot(subaccWd); accumW_dens = [accumW_dens; nanmean(accum_dens,1)]; 
    else
        subplot(subaccEd); accumE_dens = [accumE_dens; nanmean(accum_dens,1)]; 
    end 
    loglog(accum_Abin,accum_dens,'-','linewidth',1,'color',[0.75 0.75 0.75]); hold on;
    if i <=10 %ablation season
        subplot(subablWd); ablatW_dens = [ablatW_dens; nanmean(ablat_dens,1)]; 
    else
        subplot(subablEd); ablatE_dens = [ablatE_dens; nanmean(ablat_dens,1)]; 
    end 
    loglog(ablat_Abin,ablat_dens,'-','linewidth',1,'color',[0.75 0.75 0.75]); hold on;
    drawnow;
    
    %clear variables
    clear accum_Abin ablat_Abin accum_dens ablat_dens;
end
%format subplots
subplot(subablWd); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,ablatW_dens,'-','linewidth',2,'color',nanmean(ablat_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(ablatW_dens,1),'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'b) west: ablation season','FontSize',16);
ylabel('Iceberg number density'); xlabel('Iceberg surface area (m^2)'); 
subplot(subablEd); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,ablatE_dens,'-','linewidth',2,'color',nanmean(ablat_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(ablatE_dens,1),'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
set(gca,'yticklabel',[]);
text(70,0.2e-4,'d) east: ablation season','FontSize',16);
xlabel('Iceberg surface area (m^2)'); 
subplot(subaccWd); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,accumW_dens,'-','linewidth',2,'color',nanmean(accum_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(accumW_dens,1),'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'a) west: accumulation season','FontSize',16);
ylabel('Iceberg number density'); 
subplot(subaccEd); pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.15*pos(3) 1.15*pos(4)]);
loglog(F(1).Abin_mean,accumE_dens,'-','linewidth',2,'color',nanmean(accum_cmap)); hold on; grid on; 
loglog(F(1).Abin_mean,nanmean(accumE_dens,1),'-k','linewidth',2); hold on;
set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000],'xticklabel',[]); %set the x axis limits and labels
set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
text(70,0.2e-4,'c) east: accumulation season','FontSize',16);
set(gca,'yticklabel',[]);
drawnow;
saveas(modelnorm_fig,[root_dir,'Greenland-modeled-mean-seasonal-iceberg-size-density_lumped-plot.png'],'png');

% %add modeled means to observed means plots
% figure(meannorm_fig); 
% subplot(subablWn); loglog(F(1).Abin_mean,nanmean(ablatW_dens,1),'--k','linewidth',2); hold on;
% subplot(subablEn); loglog(F(1).Abin_mean,nanmean(ablatE_dens,1),'--k','linewidth',2); hold on;
% subplot(subaccWn); loglog(F(1).Abin_mean,nanmean(accumW_dens,1),'--k','linewidth',2); hold on;
% subplot(subaccEn); loglog(F(1).Abin_mean,nanmean(accumE_dens,1),'--k','linewidth',2); hold on;
% saveas(meannorm_fig,[root_dir,'Greenland-mean-seasonal-iceberg-size-density_lumped-plot.png'],'png');
clear sub*;


%% Section 4: Plot gographically-arranged observations colored by day of year
close all;

%create a HUGE composite plot in which every day of the year gets a different color to
%better resolve seasonal variations
day_sitefig = figure; set(day_sitefig,'position',[50 50 800 1200]);
%add the data to the figure template
for i = 1:length(F)
    eval(['sub',num2str(site_plotlocs(i)),' = subplot(10,2,site_plotlocs(i));']);
    
    %plot the observed size distributions
    for j = 1:size(F(i).Abin_no,1)
        loglog(F(i).Abin_mean,F(i).Abin_no(j,:),'-','linewidth',1.5,'color',day_cmap(round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))),:)); hold on;
    end
    %add labels to the annual profiles
    grid on; ylims = get(gca,'ylim');
    set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
    set(gca,'fontsize',16,'ylim',[0.1 max(ylims)],'ytick',[1e0 1e3 1e6]); %set the y axis limits
    text(100,1,[char(site_plotletts(i)),char(' '),char(F(i).site_name)],'FontSize',12);
    if site_plotlocs(i) == max(site_plotlocs)-1
%         set(gca,'xticklabel',[10^2,10^4,10^6],...
%             'yticklabel',[1 10^2 10^4 10^6]); %set the x axis limits and labels
        xlabel('Iceberg surface area (m^2)'); ylabel('Count');
    else
        set(gca,'yticklabel',[],'xticklabel',[]);
    end
    pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.25*pos(3) 1.25*pos(4)]);
    drawnow;
end
eval(['subplot(sub',num2str(max(site_plotlocs)),')']);
pos = get(gca,'position');
%plot artificial colorbar as legend
annotation('rectangle',[pos(1) pos(2)-0.06 pos(3) 0.05],'facecolor','w','edgecolor','k');
for k = 1:size(day_cmap,1)
    annotation('line',[pos(1)+17.5*(pos(3)/400)+k*(pos(3)/400) pos(1)+17.5*(pos(3)/400)+k*(pos(3)/400)],[pos(2)-0.0275 pos(2)-0.0125],'linewidth',1,'color',day_cmap(k,:)); hold on;
end
%label doy 1
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+1*(pos(3)/400)-0.01 pos(2)-0.050 0.05 0.05]; tb.String = '1'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 92
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+92*(pos(3)/400)-0.015 pos(2)-0.050 0.05 0.05]; tb.String = '92'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 183
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+183*(pos(3)/400)-0.02 pos(2)-0.050 0.05 0.05]; tb.String = '183'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 274
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+274*(pos(3)/400)-0.025 pos(2)-0.050 0.05 0.05]; tb.String = '274'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 365
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+365*(pos(3)/400)-0.03 pos(2)-0.050 0.05 0.05]; tb.String = '365'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label units
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+137*(pos(3)/400)-0.02 pos(2)-0.0635 0.15 0.05]; tb.String = 'day of year'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%add a sitemap
subplot(20,2,[site_plotlocs(site_geog_order == max(site_geog_order))+4:2:site_plotlocs(site_geog_order == max(site_geog_order)-1)-2]);
if ~exist('I')
    load('/Volumes/CALVING/Greenland_icebergs/mapping/GIMP_image_mosaic_150m.mat');
end
imagesc(I.x,I.y,I.z_adjust); axis xy equal; hold on;
colormap gray; drawnow; pos = get(gca,'position'); 
set(gca,'xlim',[-0.625e6 0.85e6],'ylim',[-3e6 -1e6],...
    'xticklabel',[],'yticklabel',[],'fontsize',16);
graticulepsn([55:5:85],[-95:10:-15]);
set(gca,'position',[0.575 0.535 1.25*pos(3) 1.25*pos(4)]);
for i = 1:length(F)
    %adjust letter plotting based on clustering & location
    if i == 1 || i == 15
        text(F(i).x+10000,F(i).y-10000,site_plotletts(i),'color','r');
    elseif i == 2 || i == 10 
        text(F(i).x+10000,F(i).y,site_plotletts(i),'color','r');
    elseif i > 10 && i ~=15
        text(F(i).x-100000,F(i).y,site_plotletts(i),'color','r');
    elseif i == 3
        text(F(i).x+10000,F(i).y,'c-d)','color','r');
    elseif i == 5
        text(F(i).x+10000,F(i).y,'e-f)','color','r');
    elseif i== 7
        text(F(i+1).x+10000,F(i+1).y,'g-i)','color','r');
    end
end
%save the figure
saveas(day_sitefig,[root_dir,'Greenland-daily-iceberg-size-distribution_subplots.png'],'png');
saveas(day_sitefig,[root_dir,'Greenland-daily-iceberg-size-distribution_subplots.eps'],'epsc');

%plot the same composite figure but only showing accumulation season data
accumday_sitefig = figure; set(accumday_sitefig,'position',[550 50 800 1200]);
%add the data to the figure template
for i = 1:length(F)
    eval(['sub',num2str(site_plotlocs(i)),' = subplot(10,2,site_plotlocs(i));']);
    
    %plot the observed size distributions
    for j = 1:size(F(i).Abin_no,1)
        if round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) <= 121 || round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) > 274
            loglog(F(i).Abin_mean,F(i).Abin_no(j,:),'-','linewidth',1.5,'color',day_cmap(round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))),:)); hold on;
        end
    end
    %add labels to the annual profiles
    grid on; ylims = get(gca,'ylim');
    set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
    set(gca,'fontsize',16,'ylim',[0.1 max(ylims)],'ytick',[1e0 1e3 1e6]); %set the y axis limits
    text(100,1,[char(site_plotletts(i)),char(' '),char(F(i).site_name)],'FontSize',12);
    if site_plotlocs(i) == max(site_plotlocs)-1
%         set(gca,'xticklabel',[10^2,10^4,10^6],...
%             'yticklabel',[1 10^2 10^4 10^6]); %set the x axis limits and labels
        xlabel('Iceberg surface area (m^2)'); ylabel('Count');
    else
        set(gca,'yticklabel',[],'xticklabel',[]);
    end
    pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.25*pos(3) 1.25*pos(4)]);
    drawnow;
end
eval(['subplot(sub',num2str(max(site_plotlocs)),')']);
pos = get(gca,'position');
%plot artificial colorbar as legend
annotation('rectangle',[pos(1) pos(2)-0.06 pos(3) 0.05],'facecolor','w','edgecolor','k');
for k = 1:size(day_cmap,1)
    annotation('line',[pos(1)+17.5*(pos(3)/400)+k*(pos(3)/400) pos(1)+17.5*(pos(3)/400)+k*(pos(3)/400)],[pos(2)-0.0275 pos(2)-0.0125],'linewidth',1,'color',day_cmap(k,:)); hold on;
end
%label doy 1
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+1*(pos(3)/400)-0.01 pos(2)-0.050 0.05 0.05]; tb.String = '1'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 92
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+92*(pos(3)/400)-0.015 pos(2)-0.050 0.05 0.05]; tb.String = '92'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 183
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+183*(pos(3)/400)-0.02 pos(2)-0.050 0.05 0.05]; tb.String = '183'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 274
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+274*(pos(3)/400)-0.025 pos(2)-0.050 0.05 0.05]; tb.String = '274'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 365
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+365*(pos(3)/400)-0.03 pos(2)-0.050 0.05 0.05]; tb.String = '365'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label units
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+137*(pos(3)/400)-0.02 pos(2)-0.0635 0.15 0.05]; tb.String = 'day of year'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%add a sitemap
subplot(20,2,[site_plotlocs(site_geog_order == max(site_geog_order))+4:2:site_plotlocs(site_geog_order == max(site_geog_order)-1)-2]);
if ~exist('I')
    load('/Volumes/CALVING/Greenland_icebergs/mapping/GIMP_image_mosaic_150m.mat');
end
imagesc(I.x,I.y,I.z_adjust); axis xy equal; hold on;
colormap gray; drawnow; pos = get(gca,'position'); 
set(gca,'xlim',[-0.625e6 0.85e6],'ylim',[-3e6 -1e6],...
    'xticklabel',[],'yticklabel',[],'fontsize',16);
graticulepsn([55:5:85],[-95:10:-15]);
set(gca,'position',[0.575 0.535 1.25*pos(3) 1.25*pos(4)]);
for i = 1:length(F)
    %adjust letter plotting based on clustering & location
    if i == 1 || i == 15
        text(F(i).x+10000,F(i).y-10000,site_plotletts(i),'color','b');
    elseif i == 2 || i == 10 
        text(F(i).x+10000,F(i).y,site_plotletts(i),'color','b');
    elseif i > 10 && i ~=15
        text(F(i).x-100000,F(i).y,site_plotletts(i),'color','b');
    elseif i == 3
        text(F(i).x+10000,F(i).y,'c-d)','color','b');
    elseif i == 5
        text(F(i).x+10000,F(i).y,'e-f)','color','b');
    elseif i== 7
        text(F(i+1).x+10000,F(i+1).y,'g-i)','color','b');
    end
end
%save the figure
saveas(accumday_sitefig,[root_dir,'Greenland-daily-winter-iceberg-size-distribution_subplots.png'],'png');


%plot the same composite figure but only showing ablation season data
ablatday_sitefig = figure; set(ablatday_sitefig,'position',[1050 50 800 1200]);
%add the data to the figure template
for i = 1:length(F)
    eval(['sub',num2str(site_plotlocs(i)),' = subplot(10,2,site_plotlocs(i));']);
    
    %plot the observed size distributions
    for j = 1:size(F(i).Abin_no,1)
        if round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) > 121 && round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) <= 274
            loglog(F(i).Abin_mean,F(i).Abin_no(j,:),'-','linewidth',1.5,'color',day_cmap(round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))),:)); hold on;
        end
    end
    %add labels to the annual profiles
    grid on; ylims = get(gca,'ylim');
    set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
    set(gca,'fontsize',16,'ylim',[0.1 max(ylims)],'ytick',[1e0 1e3 1e6]); %set the y axis limits
    text(100,1,[char(site_plotletts(i)),char(' '),char(F(i).site_name)],'FontSize',12);
    if site_plotlocs(i) == max(site_plotlocs)-1
%         set(gca,'xticklabel',[10^2,10^4,10^6],...
%             'yticklabel',[1 10^2 10^4 10^6]); %set the x axis limits and labels
        xlabel('Iceberg surface area (m^2)'); ylabel('Count');
    else
        set(gca,'yticklabel',[],'xticklabel',[]);
    end
    pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.25*pos(3) 1.25*pos(4)]);
    drawnow;
end
eval(['subplot(sub',num2str(max(site_plotlocs)),')']);
pos = get(gca,'position');
%plot artificial colorbar as legend
annotation('rectangle',[pos(1) pos(2)-0.06 pos(3) 0.05],'facecolor','w','edgecolor','k');
for k = 1:size(day_cmap,1)
    annotation('line',[pos(1)+17.5*(pos(3)/400)+k*(pos(3)/400) pos(1)+17.5*(pos(3)/400)+k*(pos(3)/400)],[pos(2)-0.0275 pos(2)-0.0125],'linewidth',1,'color',day_cmap(k,:)); hold on;
end
%label doy 1
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+1*(pos(3)/400)-0.01 pos(2)-0.050 0.05 0.05]; tb.String = '1'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 92
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+92*(pos(3)/400)-0.015 pos(2)-0.050 0.05 0.05]; tb.String = '92'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 183
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+183*(pos(3)/400)-0.02 pos(2)-0.050 0.05 0.05]; tb.String = '183'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 274
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+274*(pos(3)/400)-0.025 pos(2)-0.050 0.05 0.05]; tb.String = '274'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label doy 365
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+365*(pos(3)/400)-0.03 pos(2)-0.050 0.05 0.05]; tb.String = '365'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%label units
tb = annotation('textbox'); tb.Position = [pos(1)+17.5*(pos(3)/400)+137*(pos(3)/400)-0.02 pos(2)-0.0635 0.15 0.05]; tb.String = 'day of year'; 
tb.EdgeColor = 'none';  tb.VerticalAlignment = 'bottom'; tb.FontSize = 16;
%add a sitemap
subplot(20,2,[site_plotlocs(site_geog_order == max(site_geog_order))+4:2:site_plotlocs(site_geog_order == max(site_geog_order)-1)-2]);
if ~exist('I')
    load('/Volumes/CALVING/Greenland_icebergs/mapping/GIMP_image_mosaic_150m.mat');
end
imagesc(I.x,I.y,I.z_adjust); axis xy equal; hold on;
colormap gray; drawnow; pos = get(gca,'position'); 
set(gca,'xlim',[-0.625e6 0.85e6],'ylim',[-3e6 -1e6],...
    'xticklabel',[],'yticklabel',[],'fontsize',16);
graticulepsn([55:5:85],[-95:10:-15]);
set(gca,'position',[0.575 0.535 1.25*pos(3) 1.25*pos(4)]);
for i = 1:length(F)
    %adjust letter plotting based on clustering & location
    if i == 1 || i == 15
        text(F(i).x+10000,F(i).y-10000,site_plotletts(i),'color','r');
    elseif i == 2 || i == 10 
        text(F(i).x+10000,F(i).y,site_plotletts(i),'color','r');
    elseif i > 10 && i ~=15
        text(F(i).x-100000,F(i).y,site_plotletts(i),'color','r');
    elseif i == 3
        text(F(i).x+10000,F(i).y,'c-d)','color','r');
    elseif i == 5
        text(F(i).x+10000,F(i).y,'e-f)','color','r');
    elseif i== 7
        text(F(i+1).x+10000,F(i+1).y,'g-i)','color','r');
    end
end
%save the figure
saveas(ablatday_sitefig,[root_dir,'Greenland-daily-summer-iceberg-size-distribution_subplots.png'],'png');


%% Section 5: Plot examples for a few sites to highlight patterns
close all;

%subplots for each site: winter observed concentration, summer observed concentration, all observed & modeled density
focus_sites = [{'Upernavik South'},{'Kangerlussuaq'},{'Zachariae Isstrom'}];

%Upernavik South: highlights summer decrease in larger icebergs
%Kangerlussuaq: highlights persistent "kink"
%Zachariae: highlights tabular icebergs
example_fig = figure; set(example_fig,'position',[50 50 800 900]);
for i = 1:length(F)
    if ~isempty(strmatch(string(F(i).site_name),focus_sites))
    
    %parse the data by season
    for k = 1:length(focus_sites)
        %accumulation season iceberg concentration
        if ~isempty(strmatch(string(F(i).site_name),string(focus_sites(k))))
            accum_obsdens = []; accum_moddens = [];
            ablat_obsdens = []; ablat_moddens = [];
            
            for j = 1:size(F(i).Abin_no,1)
                %bin the data by season
                if round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) <= 121 || round(365*(F(i).dates(j,:)-floor(F(i).dates(j,:)))) > 274
                    %accumulation season ice fragment density
                    accum_obsdens = [accum_obsdens; F(i).Abin_no(j,:)./F(i).Abin_width]; 
                    accum_moddens = [accum_moddens; (F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j)))+(F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j)))];
                else
                    %ablation season iceberg concentration
                    ablat_obsdens = [ablat_obsdens; F(i).Abin_no(j,:)./F(i).Abin_width]; 
                    ablat_moddens = [ablat_moddens; (F(i).c1(j).*F(i).Abin_mean.^-F(i).alpha(j).*exp(-F(i).Abin_mean./F(i).c2(j)))+(F(i).c3(j).*exp(-F(i).Abin_mean./F(i).c4(j)))];
                end
            end
            %replace zeros with NaNs
            accum_obsdens(accum_obsdens==0) = NaN;
            ablat_obsdens(ablat_obsdens==0) = NaN;
            
            %create min & max vectors for seasonal observations
            accum_min = min(accum_obsdens,[],1); accum_max = max(accum_obsdens,[],1);
            ablat_min = min(ablat_obsdens,[],1); ablat_max = max(ablat_obsdens,[],1);
            %eliminate vector points with equal values unless for the largest iceberg size class
            match_inds = find(min(accum_obsdens,[],1) == max(accum_obsdens,[],1));
            accum_min(match_inds(1:end-1)) = NaN; clear match_inds;
            match_inds = find(min(ablat_obsdens,[],1) == max(ablat_obsdens,[],1));
            ablat_min(match_inds(1:end-1)) = NaN; clear match_inds;
                
            %plot the observed data binned by season
%             subplot(length(focus_sites),2,k*2-1);
            subplot(length(focus_sites),1,k);
            accum_lims = []; ablat_lims = [];
            accum_lims(:,1) = [F(i).Abin_mean';fliplr(F(i).Abin_mean)']; accum_lims(:,2) = [accum_min'; fliplr(accum_max)'];
            accum_lims(isnan(accum_lims(:,2)),:) = []; %remove NaNs (otherwise patch won't fill)
            p(1) = patch('Faces',[1:size(accum_lims,1)],'Vertices',accum_lims,'FaceColor',nanmean(accum_cmap),'FaceAlpha',0.5,'EdgeColor','none'); hold on;
            for j = 1:size(accum_moddens,1)
                plot(F(i).Abin_mean,accum_moddens(j,:),'-','color',nanmean(accum_cmap),'linewidth',1.5); hold on;
            end
            set(gca, 'XScale', 'log', 'YScale','log');
%             clear accum_min accum_max ablat_min ablat_max accum_lims ablat_lims p;
            
%             %create min & max vectors for seasonal models
%             accum_min = min(accum_moddens,[],1); accum_max = max(accum_moddens,[],1);
%             ablat_min = min(ablat_moddens,[],1); ablat_max = max(ablat_moddens,[],1);
%             %eliminate vector points with equal values unless for the largest iceberg size class
%             match_inds = find(min(accum_moddens,[],1) == max(accum_moddens,[],1));
%             accum_min(match_inds(1:end-1)) = NaN; clear match_inds;
%             match_inds = find(min(ablat_moddens,[],1) == max(ablat_moddens,[],1));
%             ablat_min(match_inds(1:end-1)) = NaN; clear match_inds;
            
            %plot the modeled data binned by season
%             subplot(length(focus_sites),2,k*2); %observed & modeled density
%             accum_lims = []; ablat_lims = [];
%             accum_lims(:,1) = [F(i).Abin_mean';fliplr(F(i).Abin_mean)']; accum_lims(:,2) = [accum_min'; fliplr(accum_max)'];
%             accum_lims(isnan(accum_lims(:,2)),:) = []; %remove NaNs (otherwise patch won't fill)
%             p(1) = patch('Faces',[1:size(accum_lims,1)],'Vertices',accum_lims,'FaceColor',nanmean(accum_cmap),'FaceAlpha',0.5,'EdgeColor','none'); hold on;
%             ablat_lims(:,1) = [F(i).Abin_mean';fliplr(F(i).Abin_mean)']; ablat_lims(:,2) = [ablat_min'; fliplr(ablat_max)'];
%             ablat_lims(isnan(ablat_lims(:,2)),:) = []; %remove NaNs (otherwise patch won't fill)
%             p(2) = patch('Faces',[1:size(ablat_lims,1)],'Vertices',ablat_lims,'FaceColor',nanmean(ablat_cmap),'FaceAlpha',0.5,'EdgeColor','none'); hold on;
            ablat_lims(:,1) = [F(i).Abin_mean';fliplr(F(i).Abin_mean)']; ablat_lims(:,2) = [ablat_min'; fliplr(ablat_max)'];
            ablat_lims(isnan(ablat_lims(:,2)),:) = []; %remove NaNs (otherwise patch won't fill)
            p(2) = patch('Faces',[1:size(ablat_lims,1)],'Vertices',ablat_lims,'FaceColor',nanmean(ablat_cmap),'FaceAlpha',0.3,'EdgeColor','none'); hold on;
            for j = 1:size(ablat_moddens,1)
                plot(F(i).Abin_mean,ablat_moddens(j,:),'-','color',nanmean(ablat_cmap),'linewidth',1.5); hold on;
            end
            set(gca, 'XScale', 'log', 'YScale','log');
            clear accum_min accum_max ablat_min ablat_max accum_lims ablat_lims p;
            
            %format subplots
%             subplot(length(focus_sites),2,k*2-1); grid on;
            subplot(length(focus_sites),1,k); grid on; box on; 
            set(gca,'xlim',[0 2e6],'xtick',[1e2,1e3,1e4,1e5,1e6]); %set the x axis limits and labels
            set(gca,'fontsize',16,'ylim',[1e-5 5e4],'ytick',[1e-4 1e-2 1e0 1e2 1e4]); %set the y axis limits
%             subplot(length(focus_sites),2,k*2); grid on; pos = get(gca,'position');
%             set(gca,'xlim',[0 2e6],'xtick',[100,10000,1000000]); %set the x axis limits and labels
%             set(gca,'fontsize',16,'ylim',[1e-5 2e4],'ytick',[1e-3 1e0 1e3]); %set the y axis limits
%             annotation('textbox',[pos(1) pos(2)+pos(4) pos(3) 0.025],'String',string(F(i).site_name),...
%                 'FontSize',16,'EdgeColor','none','HorizontalAlignment','center');
            text(70,0.5e-4,[char(site_plotletts(k)),' ',char(F(i).site_name)],'fontsize',16);

            %add y-axis labels
            if k == length(focus_sites)
                xlabel('Iceberg surface area (m^2)');
                ylbl = ylabel('Iceberg number density');
            else
                set(gca,'xticklabel',[]);
            end
            %add extra annotation for the first subplot
            if k == 1
            	plot([F(i).Abin_mean(2) F(i).Abin_mean(2)],[min(get(gca,'ylim')) max(get(gca,'ylim'))],'-k','linewidth',1.0); hold on; 
            end
            %adjust positions
            pos = get(gca,'position');
            set(gca,'position',[pos(1) pos(2) 1.05*pos(3) 1.25*pos(4)]);
        end
    end
    drawnow;
    
    %add additional annotation to call out important features 
    %ALL HARD-CODED IN WITH ANNOTATIONS SO THEY REQUIRE ADJUSTMENTS FOR ANY
    %CHANGES IN PLOT POSITIONS (SUBPLOT NUMBER OR SIZE)
    %tabular icebergs captured by fragmentation fits
    annotation('textbox',[0.67 0.175 0.25 0.05],'String',...
        'tabular icebergs calved via uncorrelated fracture','FontSize',16,'EdgeColor','none',...
        'HorizontalAlignment','center');
    annotation('rectangle',[0.69 0.11 0.2 0.06],'FaceColor','none','LineWidth',2);
    %potential bias in smallest "observed" icebergs
    annotation('textbox',[0.135 0.35 0.7 0.02],'String',...
        'more icebergs than predicted by fragmentation theory','FontSize',16,'EdgeColor','none');
    annotation('doublearrow',[0.13 0.13],[0.34 0.37],'HeadStyle','plain','Head1Length',6,'Head2Length',6);
    %"kink" in profiles from preferential crushing or melting
    annotation('textbox',[0.17 0.54 0.25 0.02],'String',...
        'fewer icebergs than predicted by fragmentation theory','FontSize',16,'EdgeColor','none',...
        'HorizontalAlignment','center');
    annotation('doublearrow',[0.295 0.295],[0.565 0.595],'HeadStyle','plain','Head1Length',6,'Head2Length',6);
    %seasonality in size distributions: shift to fewer larger bergs in summer
    annotation('textbox',[0.30 0.95 0.35 0.02],'String',...
        'seasonal decrease in iceberg abundance','FontSize',16,'EdgeColor','none',...
        'HorizontalAlignment','left');
    annotation('textbox',[0.15 0.95 0.15 0.02],'String',...
        'steady iceberg abundance','FontSize',16,'EdgeColor','none',...
        'HorizontalAlignment','right');
    
    end
end
%save the figure
saveas(example_fig,[root_dir,'Greenland-example-iceberg-size-distribution_subplots.png'],'png');
saveas(example_fig,[root_dir,'Greenland-example-iceberg-size-distribution_subplots.eps'],'epsc');


%% Section 6: Plot fragmentation theory-modeled size distributions
close all;

%create subplots of fit curves, with each panel corresponding to a
%different year or season
annual_fig = figure; set(annual_fig,'position',[650 50 800 1200]); 
season_fig = figure; set(season_fig,'position',[1050 50 800 1200]); 
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
        set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6],'xtick',[100,1000,10000,100000,1000000],'xticklabel',[10^2,10^3,10^4,10^5,10^6]); %set the x axis limits and labels
        set(gca,'fontsize',16,'ylim',[10^-13 1]); %set the y axis limits
        
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
        set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6],'xtick',[100,1000,10000,100000,1000000],'xticklabel',[10^2,10^3,10^4,10^5,10^6]); %set the x axis limits and labels
        set(gca,'fontsize',16,'ylim',[10^-13 1]); %set the y axis limits

        drawnow;
    end
    F(i).branching = branching; F(i).tabular = tabular; %add the expected iceberg counts to the data structure
    clear branching tabular;
end
% subplot(sub1a);
% leg = legend(pl,{num2str([2011:1:2019]')}); set(leg,'location','eastoutside');
% set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6]); xticks = get(gca,'xtick');
% set(gca,'fontsize',16,'xticklabel',xticks./10^6,'ylim',[10^-13 1]); 
% xlabel('Iceberg surface area (km^2)'); ylabel('Normalized count');
% text(0.00003*max(get(gca,'xlim')),0.4*max(get(gca,'ylim')),'a) ','fontsize',16);
% subplot(sub1b);
% legr = legend(pr,'winter','spring','summer','fall'); set(legr,'location','eastoutside');
% set(gca,'xlim',[0 ceil(max(surfA)./10^6)*10^6]); xticks = get(gca,'xtick');
% set(gca,'fontsize',16,'xticklabel',xticks./10^6,'ylim',[10^-13 1]); 
% xlabel('Iceberg surface area (km^2)'); ylabel('Normalized count');
% text(0.00003*max(get(gca,'xlim')),0.4*max(get(gca,'ylim')),'b) ','fontsize',16);

%format the annual plot figure
figure(annual_fig);
for j = 1:length([start_yr:1:end_yr])
    subplot(length([start_yr:1:end_yr]),1,j);
    text(0.00003*max(get(gca,'xlim')),0.01*max(get(gca,'ylim')),num2str(start_yr-1+j),'fontsize',16);
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
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'a) winter','fontsize',16);
subplot(4,1,2);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'b) spring','fontsize',16);
subplot(4,1,3);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'c) summer','fontsize',16);
subplot(4,1,4);
text(0.15*max(get(gca,'xlim')),0.05*max(get(gca,'ylim')),'d) fall','fontsize',16);
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
set(gca,'fontsize',16,'yticklabel',yticks./10^6); 
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
set(gca,'fontsize',16,'yticklabel',yticks./10^6); 
xlabel('Year'); ylabel('Branching size cut-off (km^2)');
clear pl pr;

%save the figures
saveas(figure2,[root_dir,'Greenland-branching-uncorrelated-calving-transition_scatterplot.eps'],'epsc');
saveas(figure2,[root_dir,'Greenland-branching-uncorrelated-calving-transition_scatterplot.png'],'png');
saveas(figure3,[root_dir,'Greenland-branching-cutoff_scatterplot.eps'],'epsc');
saveas(figure3,[root_dir,'Greenland-branching-cutoff_scatterplot.png'],'png');


