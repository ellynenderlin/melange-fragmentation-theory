%%% 3_calculate_melange_meltwater_fluxes: Use this code to combine
%%% estimated iceberg size distributions and melt rates to calculate fjord
%%% iceberg meltwater fluxes.

%% Section 0: Initialize (run every time)
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/');

%specify directories for required files ([root_dir,'/',site_abbrevs(i)])
root_dir = '/Volumes/Jokulhaup_5T/Greenland-melange/'; %include trailing / in file name
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
%% Section 1: Load iceberg size distributions & melange masks
%  TESTING FOR A SINGLE SITE... NEED TO AUTOMATE

%navigate to site folder and load data
for p = 1:length(site_abbrevs);
    cd([root_dir,site_abbrevs(p,:)]);
    load('ASG-melange-masks.mat');
    load('ASG-20110325_melange-DEMfilled.mat');
    dists = readtable('ASG-20110325-iceberg-distribution.txt');
    cd meltrates
    melts = readtable('ASG_20110325-20110411_iceberg_meltinfo.csv');
    
end

