%%% Plot figures to demonstrate application of fragmentation theory
%%% equations fit to iceberg size distributions for melange around
%%% Greenland. 

%Section 0: Initializes the code. Specify the paths for each directory and
%a few helpful variables for plotting. You must have the cmocean toolbox
%for Matlab downloaded
%(https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)
%or you can swap the call to cmocean when creating "cmap" with one of
%Matlab's built-in colormaps.

%Section 1: Uses a Landsat image and the melange-masks matlab file in each
%site directory to create an overview figure for each study site.

%Section 2: Uses example DEMs for each site that are saved in a separate
%directory, which you specify in the initialization section (Section 0),
%and the melange-masks matlab file to create example melange DEMs for each
%site. You must move example DEMs to a directory before you run this!

%% Section 0: Initialize (run every time)
clearvars; close all;

%specify directories for required files
root_dir = '/Volumes/CALVING/Greenland_icebergs/iceberg-fragmentation/'; %over-arching directory (include trailing /)
DEM_dir = root_dir; %code will navigate to this directory and then into a site name directory to look for images

%add cmocean toolbox to your Matlab path
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/');

%set plot variables
years = (2011:1:2021); %Note to self: Changed to parenthesis. 
cmap = cmocean('solar',length(years)); %color map for the terminus delineations -- Note: Changed from Haline to Parula. 
% cmap = colormap(gray(length(year))); %alternative gray-scale colormap for terminus delineations
site_names = ['HM';'KO';'AG';'IG';'UN';'US';'IB';'UM';'RI';'JI';'KB';'HH';'MG';'KL';'MD';'DJ';'ZI']; %used to identify site-specific directories (counter-clockwise from NW)
reg_flags = [3;3;3;3;3;3;2;2;2;2;1;1;1;1;4;4;5]; %specifies the region for each site listed in site_names (1=SE,2=SW,3=NW,4=CE,5=NE)
reg_colors = [215,25,28; 253,174,97; 255,255,191; 171,217,233; 44,123,182]/255; %same colorblind-friendly color scheme as analyze_iceberg_size_distribution_curve_fits.m
region_flag = ['NW';'NW';'NW';'NW';'NW';'NW';'CW';'CW';'CW';'CW';'SE';'SE';'SE';'SE';'CE';'CE';'NE']; 

%specify generic variables
rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)
WV_sigma = 2.9; %DEM uncertainty in meters (used when filtering bad data)
Wmax = 1000; %maximum likely iceberg width (based on qualitative inspection of images during terminus mapping)
Hmax = 800; %threshold thickness in meters that you do not expect icebergs to exceed (grounding line thickness is a good proxy)
cmax = 16; %melange elevation map upper limit in meters


cd(root_dir);

%% Section 1: plot each Landsat images with an overlay of the fjord polygon
disp('Plotting overview figure for each melange size distribution study site');

for i = 1:length(site_names)
    disp(site_names(i,:));
    cd([root_dir,site_names(i,:)]); %navigate to the site directory
    
    %load & plot landsat panchromatic polar stereo image
    l8 = dir('LC08*'); %find the landsat folder (file name is the same as the folder name with '_B8PS.TIF' concatenated on the end
    [I,S] = readgeoraster([l8(1).folder,'/',l8(1).name,'/',l8(1).name,'_B8PS.TIF']); %open the Landsat panchromatic image that has been reprojected to Polar Stereo coordinates
    x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
    y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
    figure; set(gcf,'position',[50 50 600 600]);
    imagesc(x,y,I); colormap gray; axis xy equal; hold on; %plot the image in grayscale
    clear l8; 
    
    %load and plot the fjord polygon
    load([site_names(i,:),'-melange-masks.mat']); %load the melange mask file
    %plot dummy lines for colorbar to ensure all years are included
    
    for j = 1:length(cmap)
       pl(j) = plot(melmask.dated(1).x,melmask.dated(1).y,'-','color',cmap(j,:),'linewidth',1.5); 
       hold on; 
       
    end
    
    %plot actual data
    for j = 1:length(melmask.dated)
       plot(melmask.dated(j).x,melmask.dated(j).y,'-','color',cmap(str2double(melmask.dated(j).datestring(1:4))-2011+1,:),'linewidth',1.5); hold on; 
    end
    %Note to self: Changed str2num to str2double 
    
    %plot the large mask for the melange and the glacier
    plot(melmask.uncropped.x,melmask.uncropped.y,'-','color','w','linewidth',2.5); hold on; %thick white line
    plot(melmask.uncropped.x,melmask.uncropped.y,'-','color',reg_colors(reg_flags(i),:),'linewidth',2); hold on; %thinner line colored to region
    
    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],...
        'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)],...
        'fontsize',20); grid on; %zoom to the extend of the bigger mask
    
    xticks = get(gca,'xtick'); yticks = get(gca,'ytick'); %get the default axis tick locations
    set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000); %change the tick labels from meters to kilometers
    xlabel('Easting (km)','fontsize',20); ylabel('Northing (km)','fontsize',20); %add axis labels
    
   

    
    %add a legend showing the colors for the terminus delineations (colors used to distinguish observation years)
    if reg_flags(i) <=2
        leg = legend(pl,num2str(years'),'location','west','orientation','vertical');
    else
        leg = legend(pl,num2str(years'),'location','east','orientation','vertical');
    end
    drawnow;
    saveas(gcf,[root_dir,site_names(i,:),'/',site_names(i,:),'_site-map.eps'],'epsc'); 
    saveas(gcf,[root_dir,site_names(i,:),'/',site_names(i,:),'_site-map.png'],'png'); %save the image
    
    hold on;
    
    disp('moving on...'); drawnow; %Took out close all so I could see them side by side, add later if needed. -A
end
disp('Move onto step 2!');


%% Section 2: plot example DEMs (only run if you have a good example DEM for each site)
disp('Plotting example DEM for each melange size distribution study site');

% %uncomment lines below to selectively replot
% clear site_names reg_flags;
% site_names = ['IG';'UN';'IB';'RI';'MD'];
% reg_flags = [3;3;2;2;4];

for i = 1:length(site_names)
    disp(site_names(i,:));
    
    %load the DEM
    cd([DEM_dir,site_names(i,:)]); %navigate to the site directory
    DEMs = dir('*dem.tif'); %identify the DEM
    DEMmat_dates = DEMs(1).name(6:13); %pull the YYYYMMDD out of the file name
    [Z.z.raw,S] = readgeoraster(DEMs(1).name); %load the DEM
    Z.x = S.XWorldLimits(1)+0.5*S.CellExtentInWorldX:S.CellExtentInWorldX:S.XWorldLimits(2)-0.5*S.CellExtentInWorldX;
    Z.y = S.YWorldLimits(2)-0.5*S.CellExtentInWorldY:-S.CellExtentInWorldY:S.YWorldLimits(1)+0.5*S.CellExtentInWorldY;
    Z.z.raw(Z.z.raw==-9999) = NaN; %replace no data value (-9999) with NaN
    
    %convert to orthometric elevations
    disp('calculating the geoid heights for the DEM...')
    warning off;
    geoid_spacing = 100; %only grab geoid height at specified increment (in meters)
    [ZXgrid,ZYgrid] = meshgrid(Z.x,Z.y); %create a grid of coordinates
    for j = 1:geoid_spacing:size(ZXgrid,1)
        for k = 1:geoid_spacing:size(ZXgrid,2)
            [Zlon(ceil(j/geoid_spacing),ceil(k/geoid_spacing)),Zlat(ceil(j/geoid_spacing),ceil(k/geoid_spacing))] = ps2wgs(ZXgrid(j,k),ZYgrid(j,k));
            G(ceil(j/geoid_spacing),ceil(k/geoid_spacing)) = geoidheight(Zlat(ceil(j/geoid_spacing),ceil(k/geoid_spacing)),Zlon(ceil(j/geoid_spacing),ceil(k/geoid_spacing)));
        end
    end
    %interpolate the sparse geoid height values to the DEM grid
    disp('converting from ellipsoidal heights to orthometric elevations...');
    Z.z.geoid = single(interp2(ZXgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),ZYgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),G,ZXgrid,ZYgrid,'linear'));
    Z.z.ortho = Z.z.raw - Z.z.geoid; %Z.z = rmfield(Z.z,'raw');

    
    %add melange masks
    load([root_dir,site_names(i,:),'/',site_names(i,:),'-melange-masks.mat']);
    disp('DEM and melange mask loaded');
    
    %convert melange mask to fjord mask
    in = inpolygon(ZXgrid,ZYgrid,melmask.uncropped.x,melmask.uncropped.y); %inpolygon identifies DEM coordinates within the large site mask
    Z.fjord.DEM_maskX = single(melmask.uncropped.x); Z.fjord.DEM_maskY = single(melmask.uncropped.y);
    Z.fjord.DEM_mask = zeros(size(Z.z.ortho)); %create a dummy boolean mask grid 
    Z.fjord.DEM_mask(in) = 1; %mask values inside the large site mask = 1
    Z.fjord.DEM_mask = round(Z.fjord.DEM_mask); Z.fjord.DEM_mask = logical(Z.fjord.DEM_mask); %convert to a logical mask
    
    %adjust sea level to zero meters to account for tides
    [Z,data_mask,gap_mask] = sl_correction(Z,WV_sigma); %sl_correction function in the same Github repo
    M.DEM.x = Z.x; M.DEM.y = Z.y;
    M.DEM.z = Z.z.adjusted; M.DEM.z(Z.z.raw==-9999) = NaN; M.DEM.z(isnan(Z.z.raw)) = NaN; 
    %isolate melange elevations
    melange = M.DEM.z; %grab all elevations from the DEM
    melange(isnan(M.DEM.z)) = 0; %replace NaNs with zeros
    melange(melange<0) = 0; melange(melange>Hmax/(917/(1026-917))) = NaN; %filter-out anomalously low & high elevations
    melange(Z.fjord.DEM_mask==0) = NaN; %replace all elevations outside the melange with NaNs
    
    %plot the melange DEM
    figDEM = figure; set(figDEM,'position',[550 100 1000 500]);
    imagesc(M.DEM.x,M.DEM.y,melange); axis xy equal; hold on; %plot the image of melange elevations
    melange_cmap = cmocean('balance',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap); %set the colormap 
    DEMax = gca; %create an axis handle
    set(gca,'clim',[0 cmax]); %set the elevation limits 
    cbar = colorbar; cbar.Label.String  = 'elevation (m a.s.l.)'; %label the colorbar
        
    
    %plot dummy lines for colorbar to ensure all years are included
    for j = 1:length(cmap)
        pl(j) = plot(melmask.dated(1).x,melmask.dated(1).y,'-','color',cmap(j,:),'linewidth',1.5); hold on; %Note to self: put zeros in front of pl(j)
    end
    %plot actual data
    for j = 1:length(melmask.dated)
        plot(melmask.dated(j).x,melmask.dated(j).y,'-','color',cmap(str2double(melmask.dated(j).datestring(1:4))-2011+1,:),'linewidth',1.5); hold on; %Note to self: Changed str2num to str2double
    end
    
    %plot the large mask for the melange and the glacier
    plot(melmask.uncropped.x,melmask.uncropped.y,'-','color','w','linewidth',2.5); hold on;
    plot(melmask.uncropped.x,melmask.uncropped.y,'-','color',reg_colors(reg_flags(i),:),'linewidth',2); hold on;
    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],...
        'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)],...
        'fontsize',20); %set the plot axis limits
    xticks = get(gca,'xtick'); yticks = get(gca,'ytick'); %get the tick labels
    set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',20); %convert tick label units from meters to kilometers
    xlabel('Easting (km)','fontsize',20); ylabel('Northing (km)','fontsize',20); %label axes
    %add a legend showing the colors for the terminus delineations (colors used to distinguish observation years)
    if reg_flags(i) <=2
        leg = legend(pl,num2str(years'),'location','west','orientation','vertical');
    else
        leg = legend(pl,num2str(years'),'location','east','orientation','vertical');
    end
    
    axpos = get(gca,'position'); set(gca,'position',[axpos(1) 0.11 axpos(3) 0.75]); axpos = get(gca,'position');  %set the position of the map in the figure window
    legpos = get(leg,'position'); %get the default position of the legend
    set(leg,'position',[(axpos(1)+axpos(3)/2)-(legpos(3)/2) (axpos(2)+axpos(4))+0.02 legpos(3) legpos(4)]); %shift the legend position
    drawnow;
    
    %save the figure
    saveas(figDEM,[root_dir,site_names(i,:),'/',site_names(i,:),'_',DEMmat_dates(:)','-melange-rawDEM.eps'],'epsc');
    saveas(figDEM,[root_dir,site_names(i,:),'/',site_names(i,:),'_',DEMmat_dates(:)','-melange-rawDEM.png'],'png');
     drawnow;
    disp('moving on...');
    clear M Z melange melmask in G Z*grid Zlat Zlon;
    
end
%% Section 3: Size distribution scatterplots
%load the aggregated data
if ~exist('F')
    load([root_dir,'Greenland-iceberg-fragmentation-curves.mat']);
end





