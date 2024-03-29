% This script runs all the functions written to perform the analysis of
% melange iceberg distributions for outlet glaciers in Greenland.

% Author:   Jukes Liu
%           Department of Geosciences
%           Boise State University
%           Date: 04/19/2021

clear all; close all;
addpath('/Users/icebergs/general-code/');
addpath('/Users/icebergs/general-code/cmocean/'); % add general functions to path)

% set paths and glacier to analyze manually:
glacier_abbrev = 'HH'; % SET GLACIER ID HERE
basepath='/Users/icebergs/iceberg-fragmentation/';
addpath([basepath,'DEMsizes_matlab-python']);
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