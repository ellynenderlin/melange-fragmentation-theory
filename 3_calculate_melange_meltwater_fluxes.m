%%% 3_calculate_melange_meltwater_fluxes: Use this code to combine
%%% estimated iceberg size distributions and melt rates to calculate fjord
%%% iceberg meltwater fluxes.

%% Section 0: Initialize (run every time)
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/LSQ_LUT_PIECEWISE');

%specify directories for required files ([root_dir,'/',site_abbrevs(i)])
% root_dir = '/Volumes/Jokulhaup_5T/Greenland-melange/'; %include trailing / in file name
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-melt/';
cd(root_dir);

%create a month naming matrix to convert 3-letter months to numeric months
month_array = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};
years = 2011:1:2020; start_yr = years(1); end_yr = years(end);

%2-letter region flagging (based on alphabetical order of site folders)
site_abbrevs = ['ASG'];
site_names = [{'Alison'}];
region_flag = ['NW'];
site_geog_order = [1]; %geographic order, counterclockwise from NW
% site_abbrevs = ['AG';'HH';'HM';'IB';'IG';'JI';'KB';'KL';'KO';'MD';'MG';'RI';'UM';'UN';'US';'ZI']; %alphabetical site directory list
% site_names = [{'Alison'},{'Helheim'},{'Ullip Sermia'},{'Salliarutsip Sermia'},{'Illulip Sermia'},...
%     {'Sermeq Kujalleq'},{'Koge Bugt'},{'Kangerlussuaq'},{'Kong Oscar'},{'Magga Dan'},...
%     {'Nigertiip Apusiia'},{'Kangilliup Sermia'},{'Umiammakku Sermia'},{'Upernavik North'},{'Upernavik South'},{'Zachariae Isstrom'}];
% region_flag = ['NW';'SE';'NW';'CW';'NW';'CW';'SE';'SE';'NW';'CE';'SE';'CW';'CW';'NW';'NW';'NE'];
% site_geog_order = [3;12;1;7;4;10;11;14;2;15;13;9;8;5;6;16]; %geographic order, counterclockwise from NW


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


disp('Initialized iceberg meltwater flux code');
%% Section 1: Convert date-specific size distribution textfiles to a concatenated csv
disp('Converting iceberg size distribution textfiles for each date to date-specific csvs & a single csv');
convert_sizedistribution_txt_to_csv(root_dir);

disp('Melange size distributions for each site saved in *-melange-distributions.csv files');
%% Section 2: Load iceberg size distributions & melange masks
close all;
%  TESTING FOR A SINGLE SITE WITH SPARSE DATA... NEED TO AUTOMATE

%navigate to site folder and load data
for p = 1:size(site_abbrevs,1);
    disp(site_abbrevs(p,:));
    cd([root_dir,site_abbrevs(p,:)]);
    load([site_abbrevs(p,:),'-melange-masks.mat']); %load all the melange DEM masks
    dists = readtable([site_abbrevs(p,:),'-melange-distributions.csv'],"VariableNamingRule","preserve"); %load all iceberg size distributions
    
    %use iceberg size distributions to estimate the submerged area of the
    %melange in each season (DJF, MAM, JJA, SON)
    dist_hdrs = dists.Properties.VariableNames;
    berg_A = dists.SurfaceArea_mean; berg_dA = dists.SurfaceArea_range; %circular surface areas (m^2)
    for j = 1:length(dist_hdrs)-2
        berg_datestr(j,:) = char(dist_hdrs(j+2));
        berg_yr(j,:) = str2num(berg_datestr(j,1:4));
        berg_mo(j,:) = str2num(berg_datestr(j,6:7));
        berg_decidate(j,:) = convert_to_decimaldate([berg_datestr(j,1:4),berg_datestr(j,6:7),berg_datestr(j,9:10)]);
    end
    berg_draft = (900/1026)*sqrt(berg_A/pi());
    %create a table of iceberg distribution numbers & normalize by area
    berg_nos = table2array(dists(:,3:end));
    berg_fracs = berg_nos./(nansum(berg_nos.*berg_A,1));
    DJF_fracs = nanmedian(berg_fracs(:,find(berg_mo==12 | berg_mo == 1 | berg_mo == 2)),2);
    MAM_fracs = nanmedian(berg_fracs(:,find(berg_mo==3 | berg_mo == 4 | berg_mo == 5)),2);
    JJA_fracs = nanmedian(berg_fracs(:,find(berg_mo==6 | berg_mo == 7 | berg_mo == 8)),2);
    SON_fracs = nanmedian(berg_fracs(:,find(berg_mo==9 | berg_mo == 10 | berg_mo == 11)),2);
    figure; set(gcf,'position',[50 50 600 600]);
    for j = 1:size(berg_fracs,2)
        pl(1) = loglog(berg_A(~isnan(berg_fracs(:,j))),berg_fracs(~isnan(berg_fracs(:,j)),j),'-','color','k'); hold on;
    end
    pl(2) = loglog(berg_A,DJF_fracs(:,1),'-b','linewidth',2); hold on;
    pl(3) = loglog(berg_A,MAM_fracs(:,1),'-g','linewidth',2); hold on;
    pl(4) = loglog(berg_A,JJA_fracs(:,1),'-m','linewidth',2); hold on;
    pl(5) = loglog(berg_A,SON_fracs(:,1),'-r','linewidth',2); hold on;
    grid on; set(gca,'fontsize',16);
    leg=legend(pl,'all','DJF','MAM','JJA','SON');
    xlabel('Iceberg surface area (m^2)','fontsize',16);
    ylabel('Count normalized by total area','fontsize',16);
    
    %loop through loading DEMs???
%     load([root_dir,site_abbrevs(p,:),'/DEMs/ASG-20110325_melange-DEMfilled.mat']);
    
    %load all meltrates & concatenate
    meltdates = []; draft = []; subA = []; dVdt = [];
    meltfiles = dir([root_dir,site_abbrevs(p,:),'/meltrates/*.csv']);
    for j = 1:length(meltfiles);
        melts = readtable([meltfiles(j).folder,'/',meltfiles(j).name]);
        draft = [draft; melts.MedianDraft_mean];  %m b.s.l.
        subA = [subA; melts.SubmergedArea_mean]; %m^2
        dVdt = [dVdt; melts.VolumeChangeRate]; %m^3/d
        meltdates = [meltdates; [repmat(string(meltfiles(j).name(length(site_abbrevs(p,:))+2:length(site_abbrevs(p,:))+9)),size(melts.MedianDraft_mean)),...
            repmat(string(meltfiles(j).name(length(site_abbrevs(p,:))+11:length(site_abbrevs(p,:))+18)),size(melts.MedianDraft_mean))]];
        clear melts;
    end
    meltrate = dVdt./subA; %m/d
    
    %visually check that the data look reasonable
    [unique_dates,refs,inds] = unique(meltdates,'rows');
    figure; set(gcf,'position',[650 50 600 600]);
    cmap = cmocean('matter',size(refs,1)); 
    for k = 1:length(inds); colors(k,:) = cmap(inds(k),:); end
    scatter(draft,meltrate,24,colors,'filled'); hold on; grid on;
    
    %fit a complex curve to the meltrate vs draft data to parameterize melt
    draft_vector = [0:1:ceil(max(berg_draft))];
    ft = fittype('poly4'); %try a 4th order polynomial to capture expected melt variability with draft
    opts = fitoptions('Method','LinearLeastSquares' );
    opts.Lower = [-Inf -Inf -Inf -Inf 0]; opts.Upper = [Inf Inf Inf Inf 0]; %force intercept through 0,0
    [fitresult,gof] = fit(draft,meltrate,ft,opts); meltcurve = feval(fitresult,draft_vector);
    plot(draft_vector(1:round(max(draft))+1),meltcurve(1:round(max(draft))+1),'-k'); hold on; %add fit curve to scatterplot of data
    %find the local maximum for deep icebergs and fix that value for all deeper icebergs
    TF = islocalmax(meltcurve(1:round(max(draft))+1)); %only find maximum over the oberved melt draft range 
    TF(round(max(draft))+2:length(draft_vector)) = 0;
    draft_meltmax = draft_vector(find(TF==1,1,'last'));
    meltrate_meltmax = meltcurve(find(TF==1,1,'last'));
    meltcurve(find(TF==1,1,'last'):end) = meltcurve(find(TF==1,1,'last'));
    plot(draft_vector(1:round(max(draft))+1),meltcurve(1:round(max(draft))+1),'--k','linewidth',2); hold on; %add adjusted fit curve to scatterplot of data
    set(gca,'fontsize',16);
    xlabel('Iceberg draft (m b.s.l.)','fontsize',16);
    ylabel('Melt rate (m/d)','fontsize',16);
    
    %estimate submerged area (read in GEEDiT melange margin delineations to get surface extents?)
    %TEMP FIX: use the median surface extent from the size distributions
    SA = nanmedian(nansum(berg_nos.*berg_A,1)); %m^2
    DJF_subA = ((2*(900/1026)+1).*berg_A).*(DJF_fracs.*SA);
    MAM_subA = ((2*(900/1026)+1).*berg_A).*(MAM_fracs.*SA);
    JJA_subA = ((2*(900/1026)+1).*berg_A).*(JJA_fracs.*SA);
    SON_subA = ((2*(900/1026)+1).*berg_A).*(SON_fracs.*SA);
    
    %multiply the surface area of each draft bin by the estimated best-fit 
    %melt rate for that draft
    berg_meltrate = interp1(draft_vector,meltcurve,berg_draft); berg_meltrate(end) = berg_meltrate(end-1);
    %binned meltwater flux (m^3/d)
    DJF_meltflux = berg_meltrate.*DJF_subA; 
    MAM_meltflux = berg_meltrate.*MAM_subA;
    JJA_meltflux = berg_meltrate.*JJA_subA;
    SON_meltflux = berg_meltrate.*SON_subA;
    
    %report seasonal meltwater fluxes in m^3/s
    fprintf('DJF typical meltwater flux = %4.0f m^3/s \n',nansum(DJF_meltflux)./86400);
    fprintf('MAM typical meltwater flux = %4.0f m^3/s \n',nansum(MAM_meltflux)./86400);
    fprintf('JJA typical meltwater flux = %4.0f m^3/s \n',nansum(JJA_meltflux)./86400);
    fprintf('SON typical meltwater flux = %4.0f m^3/s \n',nansum(SON_meltflux)./86400);
    
    %export the data to a csv table
    
    
    disp('... moving on to the next site');
end

