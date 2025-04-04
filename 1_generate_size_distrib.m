% This script runs all the functions written to perform the analysis of
% melange iceberg distributions for outlet glaciers in Greenland.

% Author:   Ellyn Enderlin & Jukes Liu
%           Department of Geosciences
%           Boise State University
%           Date: 13/09/2024

clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/');
addpath('/Users/ellynenderlin/Research/NSF_Greenland-Calving/melange-fragmentation-code/');

% set paths and glacier to analyze manually:
site_abbrev = 'DJG'; %this should be an abbreviation that is used to name your site-specific sub-directories and will become the filename prefix
basepath='/Volumes/Jokulhaup_5T/Greenland-melange/'; %this should be the overarching directory, with site-specific sub-directories
root_dir = basepath; output_dir = basepath;
cd([root_dir,site_abbrev]);

%% a) Create the melange masks using a series of manual steps
disp('Create custom melange masks for each DEM... the last step (masking) should ideally run overnight');

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

%CREATE A NEW CODE TO READ-IN A CENTERLINE (DRAW & MAKE EVEN-INCREMENT
%SPACING IN QGIS, CREATE CROSS-FJORD TRANSECTS, EXTRACT MEAN ELEVATIONS
%ALONG THE CROSS-FJORD TRANSECTS, AND DIVIDE THE MELANGE INTO MULTIPLE BINS
%THEN UPDATE THE ICEBERG DISTRIBUTION CODE TO SPIT OUT DISTRIBUTIONS FOR
%EACH BIN AND THE FULL MELANGE (CURRENTLY ONLY DOES THE FULL MELANGE)

extract_automated_iceberg_DEM_distributions(root_dir,site_abbrev,output_dir)

%% c) Automatically fit fragmentation curves to the size distributions
close all;
disp('Fit & plot size distributions... don''t close any figures while this runs');
model_size_distrib(root_dir,site_abbrev)

disp('Move on to ''2_manually_adjust_fits.ipynb'' Jupyter Notebook to correct wonky fragmentation fits');