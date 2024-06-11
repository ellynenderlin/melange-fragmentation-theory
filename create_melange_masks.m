function create_melange_masks(root_path,site_abbrev,output_path)
%------------------------------------
% Description: Intersect the melange masks used to crop digital elevation
% models with manual terminus delineations to create a time series of
% melange masks that are standardized but account for a moving boundary
%------------------------------------------

% addpath('/Users/ellynenderlin/mfiles/general/') % change this to where your general folder is

%% customization/initialization

%account for different date location in file name depending on length of site abbreviation
if length(site_abbrev) == 3
    matfile_daterefs = [5:12];
elseif length(site_abbrev) == 2
    matfile_daterefs = [4:11];
else
    error('Using a non-standard naming format! Switch to a 2- or 3-letter site abbreviation.');
end

%specify the colormap used for melange delineation
elev_cmap = cmocean('thermal',1001); elev_cmap(1,:) = [1 1 1];

%% make a mask
%if you already have a melange mask, load it, otherwise make one (requires
%a Landsat panchromatic image reprojected to Greenland Polar Stereo
%coordinates and renamed with *B8PS.TIF ending)

% find and load files
cd([output_path,'/',site_abbrev,'/']);
melmask_file = dir([site_abbrev,'*-melange-masks.mat']);
if ~isempty(melmask_file)
    disp('Loading the melange masks...');
    load([site_abbrev,'-melange-masks.mat']);
else
    disp('Create a generic fjord mask that extends over the glacier a few kilometers using a Landsat image');

    %load the panchromatic Landsat scene & create a glacier/fjord mask
    cd([root_path,'/',site_abbrev]);
    L8dir = dir('LC08*'); cd(L8dir(1).name);
    L8bands = dir('L*.TIF');
    for i = 1:length(L8bands)
        if ~isempty(strfind(L8bands(i).name,'B8PS.TIF')) %Landsat panchromatic image reprojected to PS coordinates
            if contains(version,'R2021')
                [I,S] = readgeoraster(L8bands(i).name);
            else
                [I,S] = readgeoraster(L8bands(i).name);
            end
%             info = geotiffinfo(L8bands(i).name);
            im.x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
            im.y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
            im.z = double(I);
            clear I S;
        end
    end

    %plot the figure
    figure; set(gcf,'position',[50 50 1200 1200]);
    imagesc(im.x,im.y,im.z); axis xy equal; colormap gray; drawnow; hold on;

    % plot the DEM footprints
    disp('plot DEM footprints:')
    cd([root_path,'/',site_abbrev,'/DEMs/']);
    tifs = dir('*dem.tif'); %melange DEMs provided by PGC
    for i = 1:length(tifs)
        if contains(tifs(i).name,'SETSM_')
            DEMtif_dates(i,:) = tifs(i).name(14:21); %new PGC DEM name format AFTER MOVING WV NUMBER TO NEAR THE END SO FILES SORT CHRONOLOGICALLY
        else
            DEMtif_dates(i,:) = tifs(i).name(6:13); %old PGC DEM name format
        end
    end
    DEMtif_dates = sortrows(DEMtif_dates); sizeDEMs=size(DEMtif_dates);
    melangemat_dates = '';
    melange_mats = dir([site_abbrev,'*_melange-DEM.mat']);
    for i = 1:length(melange_mats); melangemat_dates(i,:) = melange_mats(i).name(matfile_daterefs); end
    for p = 1:sizeDEMs(1)
        disp(DEMtif_dates(p,:));
        cd([root_path,'/',site_abbrev,'/']);

        %load the .tif if the .mat DEM doesn't exist
        if isempty(strmatch(DEMtif_dates(p,:),melangemat_dates))
            %load the geotiff
            if contains(version,'R2021')
                [I,S] = readgeoraster(tifs(p).name);
            else
                [I,S] = readgeoraster(tifs(p).name);
            end
            % plot the geotiff bounds
            xs = [S.XWorldLimits(1) S.XWorldLimits(1) S.XWorldLimits(2) S.XWorldLimits(2) S.XWorldLimits(1)];
            ys = [S.YWorldLimits(1) S.YWorldLimits(2) S.YWorldLimits(2) S.YWorldLimits(1) S.YWorldLimits(1)];
            plot(xs, ys,'Color','yellow'); drawnow;
        end
    end

    % take user input
    disp('Click on two vertices that define a bounding box to zoom in on the glacier and melange');
    a = ginput(2);
    set(gca,'xlim',sort(a(:,1)),'ylim',sort(a(:,2))); drawnow;

    %create the melange mask and save it
    disp('Trace fjord walls to create a fjord mask, including several km of glacier');
    [~,fjord_maskx,fjord_masky] = roipoly;
    close(gcf); drawnow;
    cd([output_path,'/',site_abbrev,'/']);
    melmask.uncropped.x = fjord_maskx; melmask.uncropped.y = fjord_masky;
    save([site_abbrev,'-melange-masks.mat'],'melmask','-v7.3');
end
melange_xlims = [min(melmask.uncropped.x)-2000 max(melmask.uncropped.x)+2000]; melange_ylims = [min(melmask.uncropped.y)-2000 max(melmask.uncropped.y)+2000];
disp('...moving on');

%% identify tifs & convert to mat-files if the DEM overlies the melange
disp('Convert geotiffs to mat-files as necessary');
melangemat_dates = ''; DEMmat_dates = ''; DEMtif_dates = '';
cd([root_path,'/',site_abbrev,'/DEMs/']);

%existing melange DEMs
melange_mats = dir([site_abbrev,'*_melange-DEM*.mat']); %include "raw" and filled versions if both were saved
for i = 1:length(melange_mats)
    melangemat_dates(i,:) = melange_mats(i).name(matfile_daterefs);
end

%full DEM matfiles from old analyses, not just melange DEM
mats = dir('*_dem.mat');
for i = 1:length(mats)
    if strcmp(site_abbrev,mats(i).name(1:2))
        DEMmat_dates(i,:) = mats(i).name(matfile_daterefs); %data saved to Y structure
    else
        DEMmat_dates(i,:) = mats(i).name(6:13); %data saved to DEM structure for crevasse analysis
    end
end

%DEM geotiffs
tifs = dir('*_dem.tif');
for i = 1:length(tifs)
    if contains(tifs(i).name,'SETSM_')
        DEMtif_dates(i,:) = tifs(i).name(14:21); %new PGC DEM name format AFTER MOVING WV NUMBER TO NEAR THE END SO FILES SORT CHRONOLOGICALLY
    else
        DEMtif_dates(i,:) = tifs(i).name(6:13); %old PGC DEM name format
    end
end


%loop through the DEM matfiles used for other research & reformat
for p = 1:length(DEMmat_dates)
    disp(['reformating DEM #',num2str(p),' of ',num2str(length(DEMmat_dates))]);
    disp(DEMmat_dates(p,:));
    cd([root_path,'/',site_abbrev,'/DEMs/']);
    DEM_name = [site_abbrev,'-',DEMmat_dates(p,:),'_melange-DEM.mat'];

    %load the .tif if the .mat DEM doesn't exist
    if isempty(strmatch(DEMmat_dates(p,:),melangemat_dates))
        %load the matfile
        load(mats(p).name);
        %account for possible alternative structure names
        if exist('Z'); Y=Z; clear Z; %from iceberg DEMs
        elseif exist('DEM'); Y=DEM; clear DEM; Y = rmfield(Y,{'xflow','yflow','zflow'}); %from crevasse DEMs
        end

        %plot a figure
        figure; set(gcf,'position',[50 50 1600 600]);
        subplot(1,2,1);
        imagesc(Y.x,Y.y,Y.z); axis xy equal; colormap(gca,elev_cmap); hold on;

        %crop the DEM
        if sign(Y.x(1)-Y.x(end)) == -1
            xmin = find(Y.x<=min(melange_xlims),1,'last'); if isempty(xmin); xmin=1; end
            xmax = find(Y.x>=max(melange_xlims),1,'first'); if isempty(xmax); xmax=length(Y.x); end
        else
            xmin = find(Y.x>=max(melange_xlims),1,'last'); if isempty(xmin); xmin=1; end
            xmax = find(Y.x<=min(melange_xlims),1,'first'); if isempty(xmax); xmax=length(Y.x); end
        end
        if sign(Y.y(1)-Y.y(end)) == -1
            ymin = find(Y.y<=min(melange_ylims),1,'last');
            ymax = find(Y.y>=max(melange_ylims),1,'first');
        else
            ymin = find(Y.y>=max(melange_ylims),1,'last');  if isempty(ymin); ymin=1; end
            ymax = find(Y.y<=min(melange_ylims),1,'first'); if isempty(ymax); ymax=length(Y.y); end
        end
        Z.x = single(Y.x(xmin:xmax)); Z.y = single(Y.y(ymin:ymax));
        Z.z.raw = single(Y.z(ymin:ymax,xmin:xmax)); Z.z.raw(Z.z.raw>3000) = NaN;
        [ZXgrid,ZYgrid] = meshgrid(Z.x,Z.y);
        subplot(1,2,2);
        imagesc(Z.x,Z.y,Z.z.raw); axis xy equal; colormap(gca,elev_cmap); hold on;
        clear Y;

        %calculate geoid heights to convert from ellipsoidal heights to orthometric elevations
        disp('calculating the geoid heights for the DEM...')
        warning off;
        geoid_spacing = 100;
        for j = 1:geoid_spacing:size(ZXgrid,1)
            for i = 1:geoid_spacing:size(ZXgrid,2)
                [Zlon(ceil(j/geoid_spacing),ceil(i/geoid_spacing)),Zlat(ceil(j/geoid_spacing),ceil(i/geoid_spacing))] = ps2wgs(ZXgrid(j,i),ZYgrid(j,i));
                G(ceil(j/geoid_spacing),ceil(i/geoid_spacing)) = geoidheight(Zlat(ceil(j/geoid_spacing),ceil(i/geoid_spacing)),Zlon(ceil(j/geoid_spacing),ceil(i/geoid_spacing)));
            end
        end
        Z.z.geoid = single(interp2(ZXgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),ZYgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),G,ZXgrid,ZYgrid,'linear'));
        disp('converting from ellipsoidal heights to orthometric elevations...');
        Z.z.ortho = Z.z.raw - Z.z.geoid; Z.z = rmfield(Z.z,'raw');
        clear Zlat Zlon G;

        %save the file if it overlaps the melange at the terminus
        close(gcf); drawnow;
        cd([root_path,'/',site_abbrev,'/DEMs/']);
        save(DEM_name,'Z','-v7.3'); %raw & intermediate elevation data
        disp(['Saved DEM geotiff from ',DEMmat_dates(p,:),' to mat-file']);
        clear Z melpoly* outline_* out_* in Z*grid;
        close all; drawnow;
    else
        disp('... already converted to correct format');
    end
end

%loop through the DEM tifs & make matfiles
for p = 1:length(tifs)
    disp(['reformating DEM tif #',num2str(p),' of ',num2str(length(tifs))]);
    disp(DEMtif_dates(p,:));
    cd([root_path,'/',site_abbrev,'/DEMs/']);
    DEM_name = [site_abbrev,'-',DEMtif_dates(p,:),'_melange-DEM.mat'];

    if isempty(strmatch(DEMtif_dates(p,:),melangemat_dates))
        %load the geotiff
        if contains(version,'R2021')
            [I,S] = readgeoraster(tifs(p).name);
        else
            [I,S] = readgeoraster(tifs(p).name);
        end
        if strcmp(S.RasterInterpretation,'cells')
            Y.x = S.XWorldLimits(1)+0.5*S.CellExtentInWorldX:S.CellExtentInWorldX:S.XWorldLimits(2)-0.5*S.CellExtentInWorldX;
            Y.y = S.YWorldLimits(2)-0.5*S.CellExtentInWorldY:-S.CellExtentInWorldY:S.YWorldLimits(1)+0.5*S.CellExtentInWorldY;
        else
            Y.x = S.XWorldLimits(1):S.CellExtentInWorldX:S.XWorldLimits(2);
            Y.y = S.YWorldLimits(2):-S.CellExtentInWorldY:S.YWorldLimits(1);
        end
        Y.z = double(I); Y.z(Y.z==-9999) = NaN;
        clear I S;

        %plot a figure
        figure; set(gcf,'position',[50 50 1600 600]);
        subplot(1,2,1);
        imagesc(Y.x,Y.y,Y.z); axis xy equal; colormap(gca,elev_cmap); hold on;

        %crop the DEM
        if sign(Y.x(1)-Y.x(end)) == -1
            xmin = find(Y.x<=min(melange_xlims),1,'last'); if isempty(xmin); xmin=1; end
            xmax = find(Y.x>=max(melange_xlims),1,'first'); if isempty(xmax); xmax=length(Y.x); end
        else
            xmin = find(Y.x>=max(melange_xlims),1,'last'); if isempty(xmin); xmin=1; end
            xmax = find(Y.x<=min(melange_xlims),1,'first'); if isempty(xmax); xmax=length(Y.x); end
        end
        if sign(Y.y(1)-Y.y(end)) == -1
            ymin = find(Y.y<=min(melange_ylims),1,'last');
            ymax = find(Y.y>=max(melange_ylims),1,'first');
        else
            ymin = find(Y.y>=max(melange_ylims),1,'last');  if isempty(ymin); ymin=1; end
            ymax = find(Y.y<=min(melange_ylims),1,'first'); if isempty(ymax); ymax=length(Y.y); end
        end
        Z.x = single(Y.x(xmin:xmax)); Z.y = single(Y.y(ymin:ymax));
        Z.z.raw = single(Y.z(ymin:ymax,xmin:xmax)); Z.z.raw(Z.z.raw>3000) = NaN;
        [ZXgrid,ZYgrid] = meshgrid(Z.x,Z.y);
        subplot(1,2,2);
        imagesc(Z.x,Z.y,Z.z.raw); axis xy equal; colormap(gca,elev_cmap); hold on;
        clear Y;

        %calculate geoid heights to convert from ellipsoidal heights to orthometric elevations
        disp('calculating the geoid heights for the DEM...')
        warning off;
        geoid_spacing = 100;
        for j = 1:geoid_spacing:size(ZXgrid,1)
            for i = 1:geoid_spacing:size(ZXgrid,2)
                [Zlon(ceil(j/geoid_spacing),ceil(i/geoid_spacing)),Zlat(ceil(j/geoid_spacing),ceil(i/geoid_spacing))] = ps2wgs(ZXgrid(j,i),ZYgrid(j,i));
                G(ceil(j/geoid_spacing),ceil(i/geoid_spacing)) = geoidheight(Zlat(ceil(j/geoid_spacing),ceil(i/geoid_spacing)),Zlon(ceil(j/geoid_spacing),ceil(i/geoid_spacing)));
            end
        end
        Z.z.geoid = single(interp2(ZXgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),ZYgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),G,ZXgrid,ZYgrid,'linear'));
        disp('converting from ellipsoidal heights to orthometric elevations...');
        Z.z.ortho = Z.z.raw - Z.z.geoid; Z.z = rmfield(Z.z,'raw');
        clear Zlat Zlon G;

        %save the file if it overlaps the melange at the terminus
        close(gcf); drawnow;
        cd([root_path,'/',site_abbrev,'/DEMs/']);
        save(DEM_name,'Z','-v7.3'); %raw & intermediate elevation data
        disp(['Saved DEM geotiff from ',DEMtif_dates(p,:),' to mat-file']);
        clear Z melpoly* outline_* out_* in Z*grid;
        close all; drawnow;
    else
        disp('... already converted to correct format');
    end
end
disp('Done with DEM matfile creation');


%% adjust the melange mask to account for coregistration differences btw Landsat & the DEMs
cd([root_path,'/',site_abbrev,'/DEMs/']);
clear melange_mats;

%re-identify DEM matfiles, focusing on only the new "raw" DEM matfiles
melangemat_dates = ''; melange_mats = dir([site_abbrev,'*_melange-DEM.mat']); %identify the new melange DEMs
for i = 1:length(melange_mats)
    melangemat_dates(i,:) = melange_mats(i).name(matfile_daterefs);
end

%edit mask as necessary
answer = questdlg('Was the fjord mask just created for the first time & has not been edited?',...
    'Mask Edit','1) Yes!','2) No!','1) Yes!');
switch answer
    case '1) Yes!'
        disp('Now modify the fjord mask using the DEMs');
        cd([root_path,'/',site_abbrev]);
        melmask_file = dir('*-melange-masks.mat');

        %load the panchromatic Landsat scene & create a glacier/fjord mask
        L8dir = dir('LC08*'); cd(L8dir(1).name);
        L8bands = dir('L*.TIF');
        for i = 1:length(L8bands)
            if ~isempty(strfind(L8bands(i).name,'B8PS.TIF')) %Landsat panchromatic image reprojected to PS coordinates
                [I,S] = readgeoraster(L8bands(i).name);
                im.x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
                im.y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
                im.z = double(I);
                clear I S;
            end
        end
        cd([root_path,'/',site_abbrev]);

        melange_xlims = [min(melmask.uncropped.x)-2000 max(melmask.uncropped.x)+2000]; melange_ylims = [min(melmask.uncropped.y)-2000 max(melmask.uncropped.y)+2000];
        %crop image to limits
        if sign(im.x(1)-im.x(end)) == -1
            xmin = find(im.x<=min(melange_xlims),1,'last'); if isempty(xmin); xmin=1; end
            xmax = find(im.x>=max(melange_xlims),1,'first'); if isempty(xmax); xmax=length(im.x); end
        else
            xmin = find(im.x>=max(melange_xlims),1,'last'); if isempty(xmin); xmin=1; end
            xmax = find(im.x<=min(melange_xlims),1,'first'); if isempty(xmax); xmax=length(im.x); end
        end
        if sign(im.y(1)-im.y(end)) == -1
            ymin = find(im.y<=min(melange_ylims),1,'last');
            ymax = find(im.y>=max(melange_ylims),1,'first');
        else
            ymin = find(im.y>=max(melange_ylims),1,'last');  if isempty(ymin); ymin=1; end
            ymax = find(im.y<=min(melange_ylims),1,'first'); if isempty(ymax); ymax=length(im.y); end
        end
        imx = single(im.x(xmin:xmax)); imy = single(im.y(ymin:ymax));
        % imz = single(im.z(ymin:ymax,xmin:xmax));
        [Xgrid,Ygrid] = meshgrid(imx,imy); %convert the vector coordinates into matrices

        %loop through the DEMs & create a DEM mosaic
        disp('looping through the DEMs & loading elevations to a 3D matrix');
        for p = 1:length(melange_mats)
            disp(['DEM #',num2str(p),' of ',num2str(length(melangemat_dates))]);
            %     disp(melangemat_dates(p,:));
            cd([root_path,'/',site_abbrev,'/DEMs/']);
            DEM_name = melange_mats(p).name;
            load(DEM_name);

            %interpolate to the cropped Landsat image
            [ZXgrid,ZYgrid] = meshgrid(Z.x,Z.y);
            DEMz(:,:,p) = interp2(ZXgrid,ZYgrid,double(Z.z.ortho),Xgrid,Ygrid); %interpolate elevations to the Landsat 8 coordinates
            clear Z*;
        end

        %plot the time-averaged DEM & redraw the melange mask
        figure1 = figure; set(gcf,'position',[50 50 1600 600]);
        imagesc(imx,imy,nanmean(DEMz,3)); axis xy equal; hold on;
        %     elev_cmap = cmocean('thermal',1001); elev_cmap(1,:) = [1 1 1];
        colormap(gca,elev_cmap); set(gca,'clim',[0 80]); cbar = colorbar;
        plot(melmask.uncropped.x,melmask.uncropped.y,'-m','linewidth',1); hold on;
        drawnow;
        disp('Redraw the fjord mask so it is custom to the DEMs');
        [~,fjord_maskx,fjord_masky] = roipoly;
        close(gcf); drawnow;

        %save the edited melange mask
        cd([output_path,'/',site_abbrev,'/']);
        melmask = rmfield(melmask,'uncropped');
        melmask.uncropped.x = fjord_maskx; melmask.uncropped.y = fjord_masky;
        save([site_abbrev,'-melange-masks.mat'],'melmask','-v7.3');
        disp('Edited mask saved');
    case '2) No!'
        disp('Carrying on with already-edited fjord mask');
end
clear answer;


%% custom modify the melange mask for each DEM (remove glacier & DEM blunders)
cd([output_path,'/',site_abbrev,'/']);
if ~exist('melmask')
    load([site_abbrev,'-melange-masks.mat']);
    disp('Fjord mask reloaded...');
end

%re-identify DEM matfiles, focusing on only the new "raw" DEM matfiles
%(this is redundant with the last subsection if you run the code straight
%through but is necessary if you have an error or otherwise break in this
%section and some new DEMs have been deleted)
cd([root_path,'/',site_abbrev,'/DEMs/']);
clear melange_mats;
melangemat_dates = ''; melange_mats = dir([site_abbrev,'*_melange-DEM.mat']); %identify the new melange DEMs
for i = 1:length(melange_mats)
    melangemat_dates(i,:) = melange_mats(i).name(matfile_daterefs);
end

%uncomment lines directly below if you are just running this section for debugging purposes
%melangemat_dates = ''; melange_mats = dir([site_abbrev,'*_melange-DEM.mat']); %identify the melange DEMs
%for i = 1:length(melange_mats)
%melangemat_dates(i,:) = melange_mats(i).name(matfile_daterefs);
%end
%tifs = dir('*_dem.tif');
%for i = 1:length(tifs)
%     if contains(tifs(i).name,'SETSM_')
%         DEMtif_dates(i,:) = tifs(i).name(14:21); %new PGC DEM name format AFTER MOVING WV NUMBER TO NEAR THE END SO FILES SORT CHRONOLOGICALLY
%     else
%         DEMtif_dates(i,:) = tifs(i).name(6:13); %old PGC DEM name format
%     end
%end

%set-up the dated melange mask structure
if ~isfield(melmask,'dated')
    melmask.dated(1).datestring = 'start'; melmask_dates(1,:) = 'empty';
else
    for j = 1:length(melmask.dated); melmask_dates(j,:) = melmask.dated(j).datestring; end
end

disp('IMPORTANT: When masking, you may need to click the hand icon on the figure just above the map & click and drag to pan!');

%trace the termini & mask-out blunders in all the DEMs
cd([root_path,'/',site_abbrev,'/DEMs/']);
for p = 1:length(melange_mats)
    %load the DEM
    disp(['DEM #',num2str(p),' of ',num2str(size(melangemat_dates,1))]);
    disp(melangemat_dates(p,:));
    DEM_name = melange_mats(p).name;
    load([root_path,'/',site_abbrev,'/DEMs/',DEM_name]);
    
    %plot the DEM
    figure1 = figure; set(gcf,'position',[50 50 1600 600]);
    imagesc(double(Z.x), double(Z.y), double(Z.z.ortho)); axis xy equal; hold on;
    colormap(gca,elev_cmap); set(gca,'clim',[0 80]); cbar = colorbar;
    plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
    title(DEM_name(1:11));
    drawnow;
    
    %identify the starting uncropped melange mask vertex to use when
    %generating cropped melange masks (should be in the water!)
    if p == 1
        [meledge_x, meledge_y] = poly2cw(melmask.uncropped.x,melmask.uncropped.y); %make sure melange outline is a clockwise polygon
        [gris_center_x, gris_center_y] = wgs2ps(-41.2, 76.7); % grab the center of the GrIS in PS coordinates
        meledge_dist = sqrt((meledge_x-gris_center_x).^2 + (meledge_y-gris_center_y).^2); %find the distance from the GrIS center to each melange outline vertex
        start_vert = find(meledge_dist==max(meledge_dist)); %find the farthest vertex & use that as the start of the outline
        set(gca,'xlim',[min([min(Z.x),min(melmask.uncropped.x)]) max([max(Z.x),max(melmask.uncropped.x)])],...
            'ylim',[min([min(Z.y),min(melmask.uncropped.y)]) max([max(Z.y),max(melmask.uncropped.y)])]);
        plot(meledge_x(start_vert(1)),meledge_y(start_vert(1)),'mx','linewidth',3); drawnow;
        
        %check that the starting vertex for the melange mask is in the ocean
        answer = questdlg('Where is the starting vertex for the melange mask (pink X)?',...
            'Mask Start Vertex','Ocean','Glacier','Ocean');
        switch answer
            case 'Ocean'
                disp('Melange mask start vertex is correct, carry on!');
            case 'Glacier'
                start_dist = sqrt((meledge_x-meledge_x(start_vert(1))).^2 + (meledge_y-meledge_y(start_vert(1))).^2);
                clear start_vert;
                start_vert = find(start_dist==max(start_dist)); %find the farthest vertex from the incorrect default starting vertex
                plot(meledge_x(start_vert(1)),meledge_y(start_vert(1)),'c+','linewidth',3); drawnow;
        end
        
        %sort melange mask vertices so that they start at the defined starting vertex in the ocean
        if start_vert(1) > 1
            outline_x = [meledge_x(start_vert(1):end); meledge_x(1:start_vert(1))];
            outline_y = [meledge_y(start_vert(1):end); meledge_y(1:start_vert(1))];
        else
            outline_x = meledge_x; outline_y = meledge_y;
        end
        clear meledge*; close(figure1); drawnow;
        
        %recreate figure without marks for the melange mask start vertex
        figure1 = figure; set(gcf,'position',[50 50 1600 600]);
        imagesc(double(Z.x),double(Z.y),double(Z.z.ortho)); axis xy equal; hold on;
        colormap(gca,elev_cmap); set(gca,'clim',[0 80]); cbar = colorbar;
        plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
        drawnow;
    end
    
    %determine if using a new DEM (if there is a geotiff with the same date
    %in the same directory)
    newtif = zeros(length(tifs),1);
    for i = 1:length(tifs)
        if ~isempty(strfind(tifs(i).name,melange_mats(p).name(matfile_daterefs)))
            newtif(i) = 1;
        end
    end
    
    %if it is a new DEM, create a custom melange mask that removes blunders, open water, and the glacier
    if sum(newtif) > 0
        clear spurious* anom* blunders;
        
        %zoom in on the DEM
        set(gca,'xlim',[max([min(melmask.uncropped.x); min(Z.x')]) min([max(melmask.uncropped.x); max(Z.x)'])],'ylim',[max([min(melmask.uncropped.y); min(Z.y')]) min([max(melmask.uncropped.y); max(Z.y')])]);
        drawnow;
        
        %check that the DEM covers a few kilometers of the melange, otherwise
        %delete the melange DEM so it isn't called again in the processing pipeline
        answer = questdlg('Does the DEM span the fjord width, cover near the terminus, & extend a few km from the terminus?',...
            'DEM coverage','1) Yes!','2) No!','1) Yes!');
        switch answer
            case '1) Yes!'
                %remind the user what to do
                disp('Iteratively generate masks to:');
                disp('   1) eliminate areas with obvious elevation blunders &')
                disp('   2) crop the melange area so it doesn''t extend over flat, open water');
                
                %iterative blunder removal
                q=1;
                while q
                    if q == 1
                        blunder_question = questdlg('Create blunder masks?',...
                            'Blunder ID','1) Yes!','2) No!','1) Yes!');
                        spurious_vals = zeros(size(Z.z.ortho));
                    else
                        blunder_question = questdlg('Are additional masks necessary?',...
                            'Blunder ID','1) Yes!','2) No!','1) Yes!');
                    end
                    %execute blunder removal based on question response
                    switch blunder_question
                        case '1) Yes!'
                            figure(figure1);
                            disp('Click on UL & LR corners of a box bounding the anomalous elevations to zoom in'); % Upper left, lower right.
                            [a] = ginput(2);
                            set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                            drawnow;
                            disp('Draw mask');
                            anom_zmask = roipoly;
                            anom_zmask = double(~anom_zmask);
                            Z.z.ortho = anom_zmask.*Z.z.ortho; Z.z.ortho(Z.z.ortho==0) = NaN;
                            imagesc(double(Z.x),double(Z.y),double(Z.z.ortho)); hold on; axis xy equal;
                            colormap(gca,elev_cmap); set(gca,'clim',[0 80]); cbar = colorbar;
                            plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
                            drawnow;
                            blunders = ~anom_zmask;
                            spurious_vals = spurious_vals + blunders;
                            q = q+1;
                            set(gca,'xlim',[max([min(melmask.uncropped.x); min(Z.x')]) min([max(melmask.uncropped.x); max(Z.x)'])],'ylim',[max([min(melmask.uncropped.y); min(Z.y')]) min([max(melmask.uncropped.y); max(Z.y')])]);
                            drawnow;
                        case '2) No!'
                            spurious_vals(spurious_vals>0) = 1;
                            clear blunders;
                            break
                    end
                    clear blunder_question;
                end
                Z.melange.blunder_mask = single(spurious_vals);
                
                %mask out the open ocean (& more blunders if discovered with the narrower elevation limits)
                q=1;
                while q
                    figure(figure1); colormap(gca,elev_cmap); set(gca,'clim',[0 9]); colorbar; drawnow;
                    if q == 1
                        mask_question = questdlg('Create open water masks?',...
                            'Open Ocean Masking','1) Yes!','2) No!','1) Yes!');
                        spurious_vals = zeros(size(Z.z.ortho));
                    else
                        mask_question = questdlg('Are additional masks necessary?',...
                            'Open Ocean Masking','1) Yes!','2) No!','1) Yes!');
                    end
                    %execute ocean masking based on question response
                    switch mask_question
                        case '1) Yes!'
                            figure(figure1);
                            disp('Click on UL & LR corners of a box bounding the open ocean to zoom in');
                            [a] = ginput(2);
                            set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                            drawnow;
                            disp('Draw mask');
                            anom_zmask = roipoly;
                            anom_zmask = double(~anom_zmask);
                            Z.z.ortho = anom_zmask.*Z.z.ortho; Z.z.ortho(Z.z.ortho==0) = NaN;
                            imagesc(double(Z.x),double(Z.y),double(Z.z.ortho)); hold on; axis xy equal;
                            colormap(gca,elev_cmap); set(gca,'clim',[0 9]); cbar = colorbar;
                            plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
                            drawnow;
                            blunders = ~anom_zmask;
                            spurious_vals = spurious_vals + blunders;
                            q = q+1;
                            set(gca,'xlim',[max([min(melmask.uncropped.x); min(Z.x')]) min([max(melmask.uncropped.x); max(Z.x)'])],'ylim',[max([min(melmask.uncropped.y); min(Z.y')]) min([max(melmask.uncropped.y); max(Z.y')])]);
                            drawnow;
                        case '2) No!'
                            spurious_vals(spurious_vals>0) = 1;
                            clear blunders;
                            break
                    end
                    clear mask_question;
                end
                Z.melange.blunder_mask = Z.melange.blunder_mask + single(spurious_vals);
                Z.melange.blunder_mask(Z.melange.blunder_mask>1) = 1;
                
                %readjust elevation color scaling
                colormap(gca,elev_cmap); set(gca,'clim',[0 80]); cbar = colorbar;
                set(gca,'xlim',[max([min(melmask.uncropped.x); min(Z.x')]) min([max(melmask.uncropped.x); max(Z.x)'])],'ylim',[max([min(melmask.uncropped.y); min(Z.y')]) min([max(melmask.uncropped.y); max(Z.y')])]);
                
                %trace the terminus, making sure to intersect the edges of the melange mask
                disp('trace a line at the bottom of the terminus cliff, placing one point on each end outside of the melange outline');
                disp('... if the terminus isn''t visible, draw a straight line across the inland limit of observations');
                [term_x,term_y,~] = improfile;
                [Z.term.x,Z.term.y] = poly2cw(term_x,term_y);
                
                %save the mask for each time step
                cd([output_path,'/',site_abbrev,'/']);
                if contains(melmask.dated(1).datestring,'start')
                    melmask.dated(1).datestring = melangemat_dates(p,:); melmask.dated(1).x = []; melmask.dated(1).y = [];
                else
                    maskref = find(contains(string(melmask_dates),melangemat_dates(p,:))==1);
                    %throw an error if any melange mask dates are
                    %duplicated in the mask structure
                    if length(maskref) > 1
                        error('Duplicate dates in the mask structure')
                    end
                    %if there is no existing date in the structure, add
                    %the date & dummy coordinate matrices to the end
                    if isempty(maskref)
                        maskref = length(melmask.dated)+1;
                    end
                    melmask.dated(maskref).datestring = melangemat_dates(p,:); melmask.dated(maskref).x = []; melmask.dated(maskref).y = [];
                end
                save([output_path,'/',site_abbrev,'/',site_abbrev,'-melange-masks.mat'],'melmask','-v7.3');
                save([root_path,'/',site_abbrev,'/DEMs/',DEM_name],'Z','-v7.3'); %raw & intermediate elevation data
            case '2) No!'
                disp('Deleting files...');
                recycle('off'); delete([root_path,'/',site_abbrev,'/DEMs/',DEM_name]); %permanently delete the DEM matfile
                DEM_sourcefiles = dir([root_path,'/',site_abbrev,'/DEMs/*',melangemat_dates(p,:),'*.t*']);
                recycle('on'); for l = 1:length(DEM_sourcefiles); delete([root_path,'/',site_abbrev,'/DEMs/',DEM_sourcefiles(l).name]); end
                clear DEM_sourcefiles;
        end
        clear answer;
        close all; drawnow;
        
    else %if the terminus was previously traced, check it to make sure it looks OK
        %zoom in on the DEM
        set(gca,'xlim',[max([min(melmask.uncropped.x); min(Z.x')]) min([max(melmask.uncropped.x); max(Z.x)'])],'ylim',[max([min(melmask.uncropped.y); min(Z.y')]) min([max(melmask.uncropped.y); max(Z.y')])]);
        drawnow;
        
        %confirm it looks good (covers the melange, doesn't look like cloud splotches in elevation
        answer = questdlg('Does the DEM extend a few km from the terminus & not look splotchy/fuzzy over the water?',...
            'DEM coverage','1) Yes!','2) No!','1) Yes!');
        switch answer
            case '1) Yes!'
                disp('Now mask bad and ocean pixels as necessary...');
                
                %iterative blunder removal
                q=1;
                while q
                    if q == 1
                        blunder_question = questdlg('Create blunder masks?',...
                            'Blunder ID','1) Yes!','2) No!','1) Yes!');
                        spurious_vals = zeros(size(Z.z.ortho));
                    else
                        blunder_question = questdlg('Are additional masks necessary?',...
                            'Blunder ID','1) Yes!','2) No!','1) Yes!');
                    end
                    %execute blunder removal based on question response
                    switch blunder_question
                        case '1) Yes!'
                            figure(figure1);
                            disp('Click on UL & LR corners of a box bounding the area to remove to zoom in');
                            [a] = ginput(2);
                            set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                            drawnow;
                            disp('Draw mask');
                            anom_zmask = roipoly;
                            anom_zmask = double(~anom_zmask);
                            Z.z.ortho = anom_zmask.*Z.z.ortho; Z.z.ortho(Z.z.ortho==0) = NaN;
                            imagesc(double(Z.x),double(Z.y),double(Z.z.ortho)); hold on; axis xy equal;
                            colormap(gca,elev_cmap); set(gca,'clim',[0 40]); cbar = colorbar;
                            plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
                            drawnow;
                            blunders = ~anom_zmask;
                            spurious_vals = spurious_vals + blunders;
                            q = q+1;
                            set(gca,'xlim',[max([min(melmask.uncropped.x); min(Z.x')]) min([max(melmask.uncropped.x); max(Z.x)'])],'ylim',[max([min(melmask.uncropped.y); min(Z.y')]) min([max(melmask.uncropped.y); max(Z.y')])]);
                            drawnow;
                        case '2) No!'
                            spurious_vals(spurious_vals>0) = 1;
                            clear blunders;
                            break
                    end
                end
                Z.melange.blunder_mask = Z.melange.blunder_mask + single(spurious_vals);
                Z.melange.blunder_mask(Z.melange.blunder_mask>1) = 1;
                
                %mask out the open ocean (& more blunders if discovered with the narrower elevation limits)
                q=1;
                while q
                    figure(figure1); colormap(gca,elev_cmap); set(gca,'clim',[0 9]); colorbar; drawnow;
                    if q == 1
                        mask_question = questdlg('Create open water masks?',...
                            'Open Ocean Masking','1) Yes!','2) No!','1) Yes!');
                        spurious_vals = zeros(size(Z.z.ortho));
                    else
                        mask_question = questdlg('Are additional masks necessary?',...
                            'Open Ocean Masking','1) Yes!','2) No!','1) Yes!');
                    end
                    %execute ocean masking based on question response
                    switch mask_question
                        case '1) Yes!'
                            figure(figure1);
                            disp('Click on UL & LR corners of a box bounding the open ocean to zoom in');
                            [a] = ginput(2);
                            set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                            drawnow;
                            disp('Draw mask');
                            anom_zmask = roipoly;
                            anom_zmask = double(~anom_zmask);
                            Z.z.ortho = anom_zmask.*Z.z.ortho; Z.z.ortho(Z.z.ortho==0) = NaN;
                            imagesc(double(Z.x),double(Z.y),double(Z.z.ortho)); hold on; axis xy equal; %%REVERT IF THIS IS WRONG - TAKE OUT DOUBLE
                            colormap(gca,elev_cmap); set(gca,'clim',[0 9]); cbar = colorbar;
                            plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
                            drawnow;
                            blunders = ~anom_zmask;
                            spurious_vals = spurious_vals + blunders;
                            q = q+1;
                            set(gca,'xlim',[max([min(melmask.uncropped.x); min(Z.x')]) min([max(melmask.uncropped.x); max(Z.x)'])],'ylim',[max([min(melmask.uncropped.y); min(Z.y')]) min([max(melmask.uncropped.y); max(Z.y')])]);
                            drawnow;
                        case '2) No!'
                            spurious_vals(spurious_vals>0) = 1;
                            clear blunders;
                            break
                    end
                    clear mask_question;
                end
                Z.melange.blunder_mask = Z.melange.blunder_mask + single(spurious_vals);
                Z.melange.blunder_mask(Z.melange.blunder_mask>1) = 1;
                
                %resave the blunder mask
                save([root_path,'/',site_abbrev,'/DEMs/',DEM_name],'Z','-v7.3'); %raw & intermediate elevation data
                
                %crop the melange mask using the terminus trace
                out_intercept = []; out_interceptx = []; out_intercepty = [];
                for i = 1:length(outline_x)-1
                    [xi,yi] = polyxpoly(outline_x(i:i+1),outline_y(i:i+1),Z.term.x,Z.term.y); %find the intersections of the terminus trace with the melange outline
                    if ~isempty(xi)
                        out_intercept = [out_intercept; i]; out_interceptx = [out_interceptx; xi]; out_intercepty = [out_intercepty; yi];
                    end
                    clear xi yi;
                end
                %find terminus positions inside the outline
                set(gca,'clim',[0 80]);
                if length(out_intercept) > 1
                    term_in = inpolygon(Z.term.x,Z.term.y,outline_x,outline_y); termx = Z.term.x(term_in); termy = Z.term.y(term_in);
                    termdist = sqrt((outline_x(out_intercept(1))-termx).^2 + (outline_y(out_intercept(1))-termy).^2);
                    term_start = find(termdist == min(termdist)); term_end = find(termdist == max(termdist));
                    if term_start > term_end
                        melpoly_x = [outline_x(1:out_intercept(1)); out_interceptx(1); termx(term_start:-1:term_end); out_interceptx(end); outline_x(out_intercept(2)+1:end)];
                        melpoly_y = [outline_y(1:out_intercept(1)); out_intercepty(1); termy(term_start:-1:term_end); out_intercepty(end); outline_y(out_intercept(2)+1:end)];
                    else
                        melpoly_x = [outline_x(1:out_intercept(1)); out_interceptx(1); termx(term_start:1:term_end); out_interceptx(end); outline_x(out_intercept(2)+1:end)];
                        melpoly_y = [outline_y(1:out_intercept(1)); out_intercepty(1); termy(term_start:1:term_end); out_intercepty(end); outline_y(out_intercept(2)+1:end)];
                    end
                    
                    %plot the cropped melange mask
                    plot(melpoly_x,melpoly_y,'--c','linewidth',2); hold on; drawnow;
                    disp('melange polygon cropped with terminus trace = cyan dashed line');
                    term_flag = 0;
                else
                    disp('old terminus trace did not intersect melange mask in 2 places');
                    plot(outline_x,outline_y,'--c','linewidth',2); hold on; drawnow;
                    term_flag = 1;
                end
                clear term_in termdist termx termy term_start term_end out_* meledge*;
                
                %retrace the terminus if it looks wonky
                disp('Terminus trace should intersect both sides of the mask and, if terminus not visible, go straight across the fjord');
                trace_question = questdlg('Does the melange mask cropped with the terminus trace look good?',...
                    'Terminus Trace Check','1) Yes!','2) No!','1) Yes!');
                switch trace_question
                    case '1) Yes!'
                        retrace_flag = 0;
                    case '2) No!'
                        retrace_flag = 1;
                end
                if retrace_flag == 1 || term_flag == 1 %retrace if manually or automatically flagged as bad
                    %trace the terminus, making sure to intersect the edges of the melange mask
                    disp('trace a line at the bottom of the terminus cliff, placing one point on each end outside of the melange outline');
                    disp('... if the terminus isn''t visible, draw a straight line across the inland limit of observations');
                    Z.term.x = []; Z.term.y = [];
                    %identify the melange mask index you should use to save the data
                    
                    maskref = find(contains(string(melmask_dates),melangemat_dates(p,:))==1);
                    if isempty(maskref)
                        if contains(melmask.dated(1).datestring,'start')
                            melmask.dated(1).datestring = melangemat_dates(p,:); melmask.dated(1).x = []; melmask.dated(1).y = [];
                        else
                            maskref = length(melmask.dated)+1;
                            melmask.dated(maskref).datestring = melangemat_dates(p,:);
                        end
                    end
                    
                    melmask.dated(maskref).x = []; melmask.dated(maskref).y = [];
                    [term_x,term_y,~] = improfile;
                    [Z.term.x,Z.term.y] = poly2cw(term_x,term_y);
                    
                    %save the mask for each time step
                    save([output_path,'/',site_abbrev,'/',site_abbrev,'-melange-masks.mat'],'melmask','-v7.3');
                    save([root_path,'/',site_abbrev,'/DEMs/',DEM_name],'Z','-v7.3'); %raw & intermediate elevation data
                    disp('new terminus trace saved');
                else
                    %identify the melange mask index you should use to save the data
                    maskref = find(contains(string(melmask_dates),melangemat_dates(p,:))==1);
                    if isempty(maskref)
                        if contains(melmask.dated(1).datestring,'start')
                            melmask.dated(1).datestring = melangemat_dates(p,:); melmask.dated(1).x = []; melmask.dated(1).y = [];
                        else
                            maskref = length(melmask.dated)+1;
                            melmask.dated(maskref).datestring = melangemat_dates(p,:);
                        end
                    end
                end
                clear prompt term_str;
                
            case '2) No!'
                disp('Deleting files...');
                recycle('off'); delete([root_path,'/',site_abbrev,'/DEMs/',DEM_name]); %permanently delete the DEM matfile
                DEM_sourcefiles = dir([root_path,'/',site_abbrev,'/DEMs/*',melangemat_dates(p,:),'*.t*']);
                recycle('on'); for l = 1:length(DEM_sourcefiles); delete([root_path,'/',site_abbrev,'/DEMs/',DEM_sourcefiles(l).name]); end
                clear DEM_sourcefiles;
        end
        clear Z melpoly* out_* in Z*grid term_flag trace_*;
        close all; drawnow;
        
    end
    
    clear newtif;
end
% clear outline_*;
disp('All anomalous elevations masked & termini delineated');
close all;

%% move all the geotiffs to a separate folder for organizational purposes (NO LONGER USED BECAUSE DEMS ARE ALREADY IN A "DEMs" FOLDER)
% 
% %create a DEMs directory as necessary
% if not(isfolder([root_path,'/',site_abbrev,'/DEMs/']))
%     mkdir([root_path,'/',site_abbrev,'/DEMs/']) % make models folder if it doesn't exist
%     disp('DEMs folder created.');
% end
% 
% %move the files
% dems = dir('*dem.tif'); orthos = dir('*ortho.tif'); metas = dir('*meta.txt');
% for p = 1:length(dems)
%    movefile(tifs(p).name,[root_path,'/',site_abbrev,'/DEMs/']);
% end
% for p = 1:length(orthos)
%    movefile(orthos(p).name,[root_path,'/',site_abbrev,'/DEMs/']);
% end
% for p = 1:length(metas)
%    movefile(metas(p).name,[root_path,'/',site_abbrev,'/DEMs/']);
% end


%% apply the fjord mask w/ moving terminus boundary
disp('Applying mask to each DEM as necessary... this takes a LONG LONG time (>1 hr/DEM)');
cd([root_path,'/',site_abbrev,'/DEMs/']);

%identify the melange DEMs (some DEMs mave have been removed in the last step due to poor coverage)
clear melangemat_dates;
melange_mats = dir([site_abbrev,'*_melange-DEM.mat']);
for i = 1:length(melange_mats)
    melangemat_dates(i,:) = melange_mats(i).name(matfile_daterefs);
end
%update the dates in the melange mask structure
clear melmask_dates;
for j = 1:length(melmask.dated); melmask_dates(j,:) = melmask.dated(j).datestring; end

% %find new DEM geotiffs (may have moved if you ran the code before and it froze part-way)
% if isfolder([root_path,'/',site_abbrev,'/DEMs/']) %if part of the code was run once & DEMs have been moved
    tifs = dir('*_dem.tif');
    for i = 1:length(tifs)
        if contains(tifs(i).name,'SETSM_')
            DEMtif_dates(i,:) = tifs(i).name(14:21); %new PGC DEM name format AFTER MOVING WV NUMBER TO NEAR THE END SO FILES SORT CHRONOLOGICALLY
        else
            DEMtif_dates(i,:) = tifs(i).name(6:13); %old PGC DEM name format
        end
    end
% else %otherwise DEMs are in root directory
%     tifs = dir('*_dem.tif');
%     for i = 1:length(tifs)
%         DEMtif_dates(i,:) = tifs(i).name(6:13);
%     end
% end

%located filled DEMs (created during next step & will be deleted if you
%decide the old mask was wrong)
filled_DEMs = dir([site_abbrev,'*melange-DEMfilled.mat']);
for i = 1:length(filled_DEMs)
    filledDEM_dates(i,:) = filled_DEMs(i).name(matfile_daterefs);
end

%loop through & mask and plot if new or just plot if old
for p = 1:length(melange_mats)
    disp(['DEM #',num2str(p),' of ',num2str(length(melangemat_dates))]);
    disp(melangemat_dates(p,:));
    cd([root_path,'/',site_abbrev,'/DEMs/']);
    DEM_name = melange_mats(p).name; load(DEM_name);

    %plot the masked DEM
    figure1 = figure; set(gcf,'position',[50 50 1600 50]);
    imagesc(double(Z.x),double(Z.y),double(Z.z.ortho)); axis xy equal;
    colormap(gca,elev_cmap); set(gca,'clim',[0 200]); cbar = colorbar; hold on;
    set(gca,'xlim',[min(Z.x) max(Z.x)],'ylim',[min(Z.y) max(Z.y)]);
    set(gcf,'position',[50 450 1600 600]);
    maskref = find(contains(string(melmask_dates),melangemat_dates(p,:))==1);
    if length(maskref) > 1
        error('Duplicate dates in the mask structure')
    else
        if ~isempty(melmask.dated(maskref).x)
            plot(melmask.dated(maskref).x,melmask.dated(maskref).y,'--c','linewidth',2); hold on;
        else
            disp('No mask was plotted because the melange mask needs to be (re)done');
        end
    end
    title(melangemat_dates(p,:),'fontsize',14); xlabel('Easting (m)','fontsize',12); ylabel('Northing (m)','fontsize',12); cbar.Label.String = 'elevation (m)';
    drawnow;

    %identify the starting uncropped melange mask vertex to use when
    %generating cropped melange masks (should be in the water!)
    if p == 1 && ~exist('start_vert')
        [meledge_x, meledge_y] = poly2cw(melmask.uncropped.x,melmask.uncropped.y); %make sure melange outline is a clockwise polygon
        [gris_center_x, gris_center_y] = wgs2ps(-41.2, 76.7); % grab the center of the GrIS in PS coordinates
        meledge_dist = sqrt((meledge_x-gris_center_x).^2 + (meledge_y-gris_center_y).^2); %find the distance from the GrIS center to each melange outline vertex
        start_vert = find(meledge_dist==max(meledge_dist)); %find the farthest vertex & use that as the start of the outline
        set(gca,'xlim',[min([min(Z.x),min(melmask.uncropped.x)]) max([max(Z.x),max(melmask.uncropped.x)])],...
            'ylim',[min([min(Z.y),min(melmask.uncropped.y)]) max([max(Z.y),max(melmask.uncropped.y)])]);
        plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
        plot(meledge_x(start_vert(1)),meledge_y(start_vert(1)),'mx','linewidth',3); drawnow;

        %check that the starting vertex for the melange mask is in the ocean
        answer = questdlg('Where is the starting vertex for the melange mask (pink X)?',...
            'Mask Start Vertex','Ocean','Glacier','Ocean');
        switch answer
            case 'Ocean'
                disp('Melange mask start vertex is correct, carry on!');
            case 'Glacier'
                start_dist = sqrt((meledge_x-meledge_x(start_vert(1))).^2 + (meledge_y-meledge_y(start_vert(1))).^2);
                clear start_vert;
                start_vert = find(start_dist==max(start_dist)); %find the farthest vertex from the incorrect default starting vertex
                plot(meledge_x(start_vert(1)),meledge_y(start_vert(1)),'c+','linewidth',3); drawnow;
        end

        %sort melange mask vertices so that they start at the defined starting vertex in the ocean
        if start_vert(1) > 1
            outline_x = [meledge_x(start_vert(1):end); meledge_x(1:start_vert(1))];
            outline_y = [meledge_y(start_vert(1):end); meledge_y(1:start_vert(1))];
        else
            outline_x = meledge_x; outline_y = meledge_y;
        end
        clear meledge*; close(figure1); drawnow;

        %recreate figure without marks for the melange mask start vertex
        figure1 = figure; set(gcf,'position',[50 50 1600 50]);
        imagesc(double(Z.x),double(Z.y),double(Z.z.ortho)); axis xy equal;
        colormap(gca,elev_cmap); set(gca,'clim',[0 200]); cbar = colorbar; hold on;
        set(gca,'xlim',[min(Z.x) max(Z.x)],'ylim',[min(Z.y) max(Z.y)]);
        set(gcf,'position',[50 450 1600 600]);
        maskref = find(contains(string(melmask_dates),melangemat_dates(p,:))==1);
        if ~isempty(melmask.dated(maskref).x)
            plot(melmask.dated(maskref).x,melmask.dated(maskref).y,'--c','linewidth',2); hold on;
        else
            disp('No mask was plotted because the melange mask needs to be redone');
        end
        title(melangemat_dates(p,:),'fontsize',14); xlabel('Easting (m)','fontsize',12); ylabel('Northing (m)','fontsize',12); cbar.Label.String = 'elevation (m)';
        drawnow;
    end

    %check if you want to forcibly re-create the mask
%     answer = questdlg('Do you need to recreate the melange mask (redrawn terminus or glacier masked)?',...
%         'Mask Redo','1) Yes!','2) No!','1) Yes!');
%     switch answer
%         case '1) Yes!'
%             redo_flag = 1;
%         case '2) No!'
%             redo_flag = 0;
%     end
%     clear answer; close(gcf);
% 
%     %identify if it is a newly-delineated terminus (if so, apply masks)
%     if isempty(melmask.dated(maskref).x) || redo_flag == 1
        %delete filled DEM created with wrong mask (if is exists)
        if ~isempty(filled_DEMs)
            %find the corresponding filled DEM
            for k = 1:size(filledDEM_dates,1)
                filledflag(k) = contains(DEM_name,filledDEM_dates(k,:));
            end
            %delete the corresponding filled DEM
            if sum(filledflag > 0)
                delete([filled_DEMs(filledflag==1).name]);
            end
            clear filledflag;
        end

        %crop the melange mask using the terminus trace
        out_intercept = []; out_interceptx = []; out_intercepty = [];
        for i = 1:length(outline_x)-1
            [xi,yi] = polyxpoly(outline_x(i:i+1),outline_y(i:i+1),Z.term.x,Z.term.y); %find the intersections of the terminus trace with the melange outline
            if ~isempty(xi)
                out_intercept = [out_intercept; i]; out_interceptx = [out_interceptx; xi]; out_intercepty = [out_intercepty; yi];
            end
            clear xi yi;
        end
        %if only one intercept redraw the terminus
        if length(out_intercept) < 2
            disp('Incorrect number of intersections between the terminus delineation w/ the fjord outline!');
            figure; set(gcf,'position',[50 50 1600 600]);
            imagesc(double(Z.x),double(Z.y),double(Z.z.ortho)); axis xy equal;
            colormap(gca,elev_cmap); set(gca,'clim',[0 80]); cbar = colorbar; hold on;
            plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
            set(gca,'xlim',[max(min(melange_xlims),min(Z.x)) min(max(melange_xlims),max(Z.x))],'ylim',[max(min(melange_ylims),min(Z.y)) min(max(melange_ylims),max(Z.y))]);
            drawnow;
            disp('retrace the bottom of the terminus cliff, placing one point on each end outside of the melange outline');
            Z.term.x = []; Z.term.y = [];
            [term_x,term_y,~] = improfile;
            [Z.term.x,Z.term.y] = poly2cw(term_x,term_y);
            close(gcf); drawnow;

            %find the intercepts with the fjord mask
            out_intercept = []; out_interceptx = []; out_intercepty = [];
            for i = 1:length(outline_x)-1
                [xi,yi] = polyxpoly(outline_x(i:i+1),outline_y(i:i+1),Z.term.x,Z.term.y); %find the intersections of the terminus trace with the melange outline
                if ~isempty(xi)
                    out_intercept = [out_intercept i]; out_interceptx = [out_interceptx xi]; out_intercepty = [out_intercepty yi];
                end
                clear xi yi;
            end
        end
        %find terminus positions inside the outline
        term_in = inpolygon(Z.term.x,Z.term.y,outline_x,outline_y); termx = Z.term.x(term_in); termy = Z.term.y(term_in);
        termdist = sqrt((outline_x(out_intercept(1))-termx).^2 + (outline_y(out_intercept(1))-termy).^2);
        term_start = find(termdist == min(termdist)); term_end = find(termdist == max(termdist));
        if term_start > term_end %start at the end of the terminus trace and move to the beginning
            melpoly_x = [outline_x(1:out_intercept(1)); out_interceptx(1); termx(end:-1:1); out_interceptx(end); outline_x(out_intercept(end)+1:end)];
            melpoly_y = [outline_y(1:out_intercept(1)); out_intercepty(1); termy(end:-1:1); out_intercepty(end); outline_y(out_intercept(end)+1:end)];
        else %start at the beginning of the terminus trace and move to the end
            melpoly_x = [outline_x(1:out_intercept(1)); out_interceptx(1); termx(1:1:end); out_interceptx(end); outline_x(out_intercept(end)+1:end)];
            melpoly_y = [outline_y(1:out_intercept(1)); out_intercepty(1); termy(1:1:end); out_intercepty(end); outline_y(out_intercept(end)+1:end)];
        end
        clear term_in termdist termx termy term_start term_end;
        figure; set(gcf,'position',[50 50 1600 600]);
        imagesc(double(Z.x),double(Z.y),double(Z.z.ortho)); axis xy equal;
        colormap(gca,elev_cmap); set(gca,'clim',[0 200]); cbar = colorbar; hold on;
        set(gca,'xlim',[max(min(melpoly_x),min(Z.x)) min(max(melpoly_x),max(Z.x))],'ylim',[max(min(melpoly_y),min(Z.y)) min(max(melpoly_y),max(Z.y))]);
        plot(melpoly_x,melpoly_y,'--c','linewidth',2); hold on;
        title(melangemat_dates(p,:),'fontsize',14); xlabel('Easting (m)','fontsize',12); ylabel('Northing (m)','fontsize',12); cbar.Label.String = 'elevation (m)';
        drawnow;

        %create melange mask
        disp('Creating a BW mask from the melange mask shapefile');
        [ZXgrid,ZYgrid] = meshgrid(Z.x,Z.y);
        in = inpolygon(ZXgrid,ZYgrid,melpoly_x,melpoly_y);
        Z.fjord.DEM_maskX = single(melpoly_x); Z.fjord.DEM_maskY = single(melpoly_y);
        Z.fjord.DEM_mask = zeros(size(Z.z.ortho));
        Z.fjord.DEM_mask(in) = 1;
        Z.fjord.DEM_mask = round(Z.fjord.DEM_mask); Z.fjord.DEM_mask = logical(Z.fjord.DEM_mask);
        set(gca,'xlim',sort([Z.x(find(sum(Z.fjord.DEM_mask,1)>=1,1,'first')) Z.x(find(sum(Z.fjord.DEM_mask,1)>=1,1,'last'))]),'ylim',sort([Z.y(find(sum(Z.fjord.DEM_mask,2)>=1,1,'first')) Z.y(find(sum(Z.fjord.DEM_mask,2)>=1,1,'last'))]));
        set(gca,'clim',[0 70]); grid on;

        %save the mask for each time step
        melmask.dated(maskref).datestring = melangemat_dates(p,:);
        melmask.dated(maskref).x = melpoly_x; melmask.dated(maskref).y = melpoly_y;
        cd([output_path,'/',site_abbrev,'/']);
        save([site_abbrev,'-melange-masks.mat'],'melmask','-v7.3');
        save([root_path,'/',site_abbrev,'/DEMs/',DEM_name],'Z','-v7.3'); %raw & intermediate elevation data
        saveas(gcf,[site_abbrev,'-',melangemat_dates(p,:),'-melange-DEMmap.png'],'png');
        disp(['Saved ',DEM_name]);
        clear Z melpoly* out_* *in Z*grid;
%     end
    close all; drawnow; clear newtif;
end

%show all the melange outlines cropped to the terminus as a quality check
%(use fix_individual_melange_masks.m as needed)
figure;
plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',2); axis xy equal; hold on;
melmask_cmap = colormap(cool(length(melmask.dated)));
for p = 1:length(melmask.dated)
    plot(melmask.dated(p).x,melmask.dated(p).y,'-','color',melmask_cmap(p,:)); hold on;
end
leg = legend; set(leg,'location','eastoutside','NumColumns',2);
drawnow;

disp('Done converting geotiffs to mat-files... move on to next function!');



end
