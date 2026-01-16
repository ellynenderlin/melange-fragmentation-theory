%% Create gifs of satellite images
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/');
% download S2 images with preautorift then make gifs
%identfy the files (run this section for each study site, specifying
%site-specific parameters below before each rerun)

%site-specific info
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';
site_abbrev = 'ASG'; site_name = 'Alison';
ref_image = 'S2A_21XWC_20180304_0_L2A_B08_clipped.tif'; %Alison = 'S2A_21XWC_20200731_1_L2A_B08_clipped.tif', %Zachariae = 'S2A_27XWH_20200802_3_L2A_B08_clipped.tif'

%load the data
load([root_dir,site_abbrev,'/',site_abbrev,'-melange-masks.mat']); %load the melange mask file
ims = dir([root_dir,site_abbrev,'/images/S2/','S*B08_clipped.tif']); im_refs = []; im_dates = [];
for k = 1:length(ims)
    % if contains(ims(k).name,'_2019') || contains(ims(k).name,'_2020')
    if contains(ims(k).name,'_2018')
%         ref_image = [ims(k).folder,'/',ims(k).name];
        im_refs = [im_refs; k]; im_dates = [im_dates;ims(k).name(11:18);];
    end
end

%sort the images by date
for k = 1:length(im_refs)
    decidate(k,1) = convert_to_decimaldate(im_dates(k,:)); 
end
[sorted_dates,date_refs] = sort(decidate);

%load a good reference image from late July 2020 (time period used as the
%reference for terminus delineations)
[I,R] = readgeoraster([root_dir,site_abbrev,'/images/S2/',ref_image]);
im.x = linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(2));
im.y = linspace(R.YWorldLimits(2),R.YWorldLimits(1),R.RasterSize(1));
im.z = double(I);
clear I R;

%create the map template with the reference date shown to start
map_fig = figure; set(map_fig,'position',[850 50 800 600]);
imagesc(im.x,im.y,imadjust(im.z./max(max(im.z)))); axis xy equal; colormap gray; drawnow; hold on;
xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
set(gca,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
title(site_name);
clear im;

%loop through the images and create a gif
nimages = 1;
for j = 1:length(date_refs)
    %decide whether to plot the image
    if j == 1
        %load the image for the specified date
        [I,R] = readgeoraster([root_dir,site_abbrev,'/images/S2/',ims(im_refs(date_refs(j))).name]);
        im.x = linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(2));
        im.y = linspace(R.YWorldLimits(2),R.YWorldLimits(1),R.RasterSize(1));
        im.z = double(I);
        clear I R;
        %plot
        imagesc(im.x,im.y,imadjust(im.z./max(max(im.z)))); axis xy equal; colormap gray; drawnow; hold on;
        set(gca,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
        xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        title([im_dates(im_refs(date_refs(j)),1:4),'/',im_dates(im_refs(date_refs(j)),5:6),'/',im_dates(im_refs(date_refs(j)),7:8)]);
        drawnow;
        last_date = sorted_dates(j);
        
        frame = getframe(map_fig);
        gif_im{nimages} = frame2im(frame); nimages = nimages+1;
    else
        if sorted_dates(j)-last_date >= 7/365
            %load the image for the specified date
            [I,R] = readgeoraster([root_dir,site_abbrev,'/images/S2/',ims(im_refs(date_refs(j))).name]);
            im.x = linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(2));
            im.y = linspace(R.YWorldLimits(2),R.YWorldLimits(1),R.RasterSize(1));
            im.z = double(I);
            clear I R;
            %plot
            imagesc(im.x,im.y,imadjust(im.z./max(max(im.z)))); axis xy equal; colormap gray; drawnow; hold on;
            set(gca,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
            xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
            title([im_dates(im_refs(date_refs(j)),1:4),'/',im_dates(im_refs(date_refs(j)),5:6),'/',im_dates(im_refs(date_refs(j)),7:8)]);
            drawnow;
            last_date = sorted_dates(j);
            
            frame = getframe(map_fig);
            gif_im{nimages} = frame2im(frame); nimages = nimages+1;
        end
    end
    
    %     cla;
    %     title('');
    %     frame = getframe(map_fig);
    %     gif_im{nimages} = frame2im(frame); nimages = nimages+1;
    clear im;
end
close;

filename = [root_dir,site_abbrev,'/',site_abbrev,'-Sentinel2-images.gif']; % Specify the output file name
for idx = 1:nimages-1
    [A,map] = rgb2ind(gif_im{idx},256);
    if idx == 1
        imwrite(A,map,filename,"gif",LoopCount=Inf, ...
            DelayTime=1)
    else
        imwrite(A,map,filename,"gif",WriteMode="append", ...
            DelayTime=1)
    end
end
close all; clear A frame gif_im map;
clear sort* im_* ims date_refs decidate;