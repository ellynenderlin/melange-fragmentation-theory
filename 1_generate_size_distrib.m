%%% This script runs all the functions written to perform the analysis of
%%% melange iceberg distributions for outlet glaciers in Greenland.

% Author:   Ellyn Enderlin & Jukes Liu
%           Department of Geosciences
%           Boise State University
%           Date: 27/05/2025

clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/inpoly2/');
addpath('/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-fragmentation-code/');

% set paths and glacier to analyze manually:
site_abbrev = 'KOG'; %this should be an abbreviation that is used to name your site-specific sub-directories and will become the filename prefix
basepath='/Volumes/Jokulhaup_5T/Greenland-melange/'; %this should be the overarching directory, with site-specific sub-directories
root_dir = basepath; output_dir = [root_dir,site_abbrev,'/'];
LCdir = dir([root_dir,site_abbrev,'/LC*']); im_dir = [LCdir(1).folder,'/',LCdir(1).name,'/']; %Landsat 8 or 9 unzipped image directory for mapping
vel_dir = '/Users/ellynenderlin/Research/miscellaneous/Greenland-VelMosaic_1995-2015/'; %GrIMP velocity mosaic: https://nsidc.org/grimp
cd([root_dir,site_abbrev]);
disp('Paths set, move along!');

%% a) Create the melange masks using a series of manual steps
disp('Create custom melange masks for each DEM... the last step (masking) might take >1 hr if there are multiple DEMs');

%check DEM location
answer = questdlg('Are all DEM matfiles & tifs located in a DEMs directory?',...
    'DEM check','1) Yes','2) No','1) Yes');
switch answer
    case '1) Yes'
        disp('carry on!');
    case '2) No'
        disp('Move files then enter "dbcont" in the command line to continue');
        keyboard
end
%create masks
create_melange_masks(root_dir,site_abbrev,output_dir) 


disp('Now you can delete the DEM geotiffs to save space on your computer!');
%% b) Extract the automated iceberg distributions from the DEM
disp('Extract elevation profiles & iceberg size distributions...');

%load the melange masks if you skipped over the last section because masks
%were previously created
w = who;
if sum(contains(w,'melmask')) == 0; load([output_dir,'/',site_abbrev,'-melange-masks.mat']); end

%manually delineate the centerline: requires a reference Landsat image & velocity vector (vx,vy) geotiffs
answer = questdlg('Do centerline & cross-flow profiles exist?',...
    'profile creation','1) Yes','2) No','2) No');
switch answer
    case '1) Yes'
        if exist('AF') == 1
            disp('profiles in workspace, advancing to data extraction...');
        else
            disp('loading existing shapefiles');
            warning('off', 'map:shapefile:missingDBF');
            shp_files = dir([root_dir,site_abbrev,'/shapefiles/',site_abbrev,'*.shp']);
            for j = 1:length(shp_files)
                if contains(shp_files(j).name,'centerline')
                    AF = shaperead([shp_files(j).folder,'/',shp_files(j).name]);
                elseif contains(shp_files(j).name,'transects')
                    XF = shaperead([shp_files(j).folder,'/',shp_files(j).name]);
                end
            end
        end
    case '2) No'
        disp('creating a centerline profile & evenly-spaced cross-flow transects...');
        [AF,XF] = create_profile_and_transects([root_dir,site_abbrev,'/'],melmask,im_dir,vel_dir,3413);
end
clear answer;

%automatically extract size distributions for the full melange and subsets
%that are defined by the transects
extract_automated_iceberg_DEM_distributions(root_dir,site_abbrev,AF,XF,im_dir,output_dir)

%% c) Automatically fit fragmentation curves to the size distributions
close all;
disp('Fit & plot size distributions... don''t close any figures while this runs');
model_size_distrib(root_dir,site_abbrev)

disp('Move on to ''2_manually_adjust_fits.ipynb'' Jupyter Notebook to correct wonky fragmentation fits');