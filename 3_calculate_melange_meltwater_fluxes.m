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
%  TESTING FOR A SINGLE SITE WITH SPARSE DATA... NEED TO AUTOMATE

%navigate to site folder and load data
for p = 1:size(site_abbrevs,1);
    cd([root_dir,site_abbrevs(p,:)]);
    load([site_abbrevs(p,:),'-melange-masks.mat']); %load all the melange DEM masks
    dists = readtable([site_abbrevs(p,:),'-melange-distributions.csv'],"VariableNamingRule","preserve"); %load all iceberg size distributions
    
    %use iceberg size distributions to estimate the submerged area of the
    %melange in each season (DJF, MAM, JJA, SON)
    dist_hdrs = dists.Properties.VariableNames;
    berg_A = dists.SurfaceArea_mean; berg_dA = dists.SurfaceArea_range;
    for j = 1:length(dist_hdrs)-2
        berg_datestr(j,:) = char(dist_hdrs(j+2));
        berg_yr(j,:) = str2num(berg_datestr(j,1:4));
        berg_mo(j,:) = str2num(berg_datestr(j,6:7));
        berg_decidate(j,:) = convert_to_decimaldate([berg_datestr(j,1:4),berg_datestr(j,6:7),berg_datestr(j,9:10)]);
    end
    DJF_nos = nanmedian(table2array(dists(:,find(berg_mo==12 | berg_mo == 1 | berg_mo == 2)+2)),2);
    MAM_nos = nanmedian(table2array(dists(:,find(berg_mo==3 | berg_mo == 4 | berg_mo == 5)+2)),2);
    JJA_nos = nanmedian(table2array(dists(:,find(berg_mo==6 | berg_mo == 7 | berg_mo == 8)+2)),2);
    SON_nos = nanmedian(table2array(dists(:,find(berg_mo==9 | berg_mo == 10 | berg_mo == 11)+2)),2);
    
    
    %loop through loading DEMs???
%     load([root_dir,site_abbrevs(p,:),'/DEMs/ASG-20110325_melange-DEMfilled.mat']);
    
    %load all meltrates & concatenate
    meltdates = []; draft = []; subA = []; dVdt = [];
    meltfiles = dir([root_dir,site_abbrevs(p,:),'/meltrates/updated/*.csv']);
    for j = 1:length(meltfiles);
        melts = readtable([meltfiles(j).folder,'/',meltfiles(j).name]);
        draft = [draft; melts.MedianDraft_mean]; 
        subA = [subA; melts.SubmergedArea_mean]; 
        dVdt = [dVdt; melts.VolumeChangeRate];
        meltdates = [meltdates; [repmat(string(meltfiles(j).name(length(site_abbrevs(p,:))+2:length(site_abbrevs(p,:))+9)),size(melts.MedianDraft_mean)),...
            repmat(string(meltfiles(j).name(length(site_abbrevs(p,:))+11:length(site_abbrevs(p,:))+18)),size(melts.MedianDraft_mean))]];
        clear melts;
    end
    meltrate = dVdt./subA;
    
    %visually check that the data look reasonable
    [unique_dates,refs,inds] = unique(meltdates,'rows');
    figure; cmap = cmocean('matter',size(refs,1)); 
    for k = 1:length(inds); colors(k,:) = cmap(inds(k),:); end
    scatter(draft,meltrate,24,colors,'filled'); hold on; grid on;
    
    %fit a complex curve to the meltrate vs draft data to parameterize melt
    draft_vector = [0:1:max(draft)];
    ft = fittype('poly4'); %try a 4th order polynomial to capture expected melt variability with draft
    opts = fitoptions('Method','LinearLeastSquares' );
    opts.Lower = [-Inf -Inf -Inf -Inf 0]; opts.Upper = [Inf Inf Inf Inf 0]; %force intercept through 0,0
    [fitresult,gof] = fit(draft,meltrate,ft,opts); meltcurve = feval(fitresult,draft_vector);
    plot(draft_vector,meltcurve,'-k'); hold on; %add fit curve to scatterplot of data
    %find the local maximum for deep icebergs and fix that value for all deeper icebergs
    TF = islocalmax(meltcurve); 
    draft_meltmax = draft_vector(find(TF==1,1,'last'));
    meltrate_meltmax = meltcurve(find(TF==1,1,'last'));
    meltcurve(find(TF==1,1,'last'):end) = meltcurve(find(TF==1,1,'last'));
    plot(draft_vector,meltcurve,'--k','linewidth',2); hold on; %add adjusted fit curve to scatterplot of data
    
    
    
end

