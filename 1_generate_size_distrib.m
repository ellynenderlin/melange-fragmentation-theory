% This script runs all the functions written to perform the analysis of
% melange iceberg distributions for outlet glaciers in Greenland.

% Author:   Jukes Liu
%           Department of Geosciences
%           Boise State University
%           Date: 04/19/2021

clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/');
addpath('/Users/ellynenderlin/Research/NSF_Greenland-Calving/iceberg-calving/DEMsizes_matlab-python');

% set paths and glacier to analyze manually:
glacier_abbrev = 'ASG'; %this should be an abbreviation that is used to name your site-specific sub-directories and will become the filename prefix
basepath='/Volumes/Jokulhaup_5T/Greenland-melange/'; %this should be the overarching directory, with site-specific sub-directories
root_path = basepath; output_path = basepath;
cd([root_path,glacier_abbrev]);

%% a) Create the melange masks using a series of manual steps
disp('Create custom melange masks... the last step (masking) should ideally run overnight');
create_melange_masks(root_path,glacier_abbrev,output_path) 

%% b) Extract the automated iceberg distributions from the DEM
extract_automated_iceberg_DEM_distributions(root_path,glacier_abbrev,output_path)

%% c) Automatically fit fragmentation curves to the size distributions
close all;
disp('Fit & plot size distributions... don''t close any figures while this runs');
model_size_distrib(root_path,glacier_abbrev)

disp('Move on to ''2_manually_adjust_fits.ipynb'' Jupyter Notebook to correct wonky fragmentation fits');