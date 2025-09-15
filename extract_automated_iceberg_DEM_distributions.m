function extract_automated_iceberg_DEM_distributions(root_dir,site_abbrev,AF,XF,mask_check,transect_spacer,im_dir,output_dir)
%------------------------------------
% Description: Construct iceberg size distributions from
% automatically-extracted elevation distributions. Saves size distributions
% to the site-specific "melange-DEMfilled.mat" file and to a tab-delimited
% text file ([site_abbrev,'-',num2str(DEM_dates(p,:)),'-iceberg-distribution.txt']).
%
% Inputs:
% root_dir = overarching directory containing site-specific directories
% site_abbrev = 2- or 3-letter site abbreviation (should be the directory
% name and is used for output file names)
% AF = along-flow profile created from create_profile_and_transects.m
% XF = across-flow transects created from create_profile_and_transects.m
% mask_check = automatically filter-out small polygons that were skipped in
% manual masking if mask_check = 1
% output_dir = directory where data products produced from the pipeline are
% stored
%
% may need to add some paths (depending on permissions) before running: 
% run in the command window when specifying root_dir, site_abbrev,output_dir
% addpath('/users/ellynenderlin/mfiles','/users/ellynenderlin/mfiles/general',...
%     '/users/ellynenderlin/mfiles/Greenland-calving/icebergs',...
%     '/users/ellynenderlin/mfiles/Greenland-calving/icebergs/automated-iceberg-mapping',...
%     '/users/ellynenderlin/mfiles/Greenland-calving/icebergs/fragtheory_data_fit');
%
%------------------------------------------

%% INITIALIZE

%account for different date location in file name depending on length of site abbreviation
if length(site_abbrev) == 3
    matfile_daterefs = [5:12];
elseif length(site_abbrev) == 2
    matfile_daterefs = [4:11];
else
    error('Using a non-standard naming format! Switch to a 2- or 3-letter site abbreviation.');
end

% specify variables
rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)
WV_sigma = 2.9; %DEM uncertainty in meters (used when filtering bad data)
z_inc = ceil(WV_sigma); %increment for melange elevation bins (min ~= WV_sigma; controls number of iceberg size classes)
threshold_range = [z_inc:1:12]; %test range for elevation thresholds to delineate icebergs using elevation contours from automated_manual_iceberg_aspects_comparison.m (needed for plot colors)
adjuster = z_inc-1; %scalar used to assign outputs of thresholding to a structure starting at index 1
Wmax = 1000; %maximum likely iceberg width (based on qualitative inspection of images during terminus mapping)
Hmax = 800; %threshold thickness in meters that you do not expect icebergs to exceed (grounding line thickness is a good proxy)
ARcomp.best.autoALL = 2; % iceberg aspect ratio 
ARcomp.range.autoALL = [1.7, 2.3]; % range in iceberg aspect ratio

% if you want to overwrite the data in the filled DEMs that were
% previously created, switch filled_flag to 1
filled_flag = 0; 

%load data
cd(output_dir);
load([site_abbrev,'-melange-masks.mat']); %created using create_melange_masks
mask_dates = '';
for i = 1:length(melmask.dated)
    mask_dates(i,:) = melmask.dated(i).datestring;
end

%extract decimal dates for DEMs
cd([root_dir,'/',site_abbrev,'/DEMs/']);
mats = dir([site_abbrev,'*melange-DEM.mat']);
leap_doys = [31 29 31 30 31 30 31 31 30 31 30 31]; norm_doys = [31 28 31 30 31 30 31 31 30 31 30 31];
leap_cumdoys = cumsum(leap_doys); leap_cumdoys = [0 leap_cumdoys(1:11)];
norm_cumdoys = cumsum(norm_doys); norm_cumdoys = [0 norm_cumdoys(1:11)];
for i = 1:length(mats)
    DEMmat_dates(i,:) = mats(i).name(matfile_daterefs);
    if mod(str2num(DEMmat_dates(1:4)),4) ~= 0
        DEM_decidate(i,:) = str2num(DEMmat_dates(i,1:4)) + (norm_cumdoys(str2num(DEMmat_dates(i,5:6)))+str2num(DEMmat_dates(i,7:8)))./sum(norm_doys);
    else
        DEM_decidate(i,:) = str2num(DEMmat_dates(i,1:4)) + (leap_cumdoys(str2num(DEMmat_dates(i,5:6)))+str2num(DEMmat_dates(i,7:8)))./sum(leap_doys);
    end
end
filled_DEMs = dir([site_abbrev,'*melange-DEMfilled.mat']);
if isempty(filled_DEMs)
    filledDEM_dates = ''; %empty string if no filled DEMs exist
else
    for i = 1:length(filled_DEMs)
        filledDEM_dates(i,:) = filled_DEMs(i).name(matfile_daterefs);
    end
end

%identify existing melange iceberg distribution datasets
cd(output_dir);
% melange_mats = dir([site_abbrev,'*_melange-distribution.mat']); melange_dates = ''; %manual iceberg sizes
% for i = 1:length(melange_mats)
%     melange_dates(i,:) = melange_mats(i).name(matfile_daterefs);
% end
iceberg_mats = dir([site_abbrev,'*_iceberg-data.mat']); iceberg_dates = ''; %automated iceberg sizes
for i = 1:length(iceberg_mats)
    iceberg_dates(i,:) = iceberg_mats(i).name(matfile_daterefs);
end

alphabet = 'abcdefghijklmnopqrstuvwxyz';

disp('Initialized size distribution extraction process.');
%% PROCESS MELANGE MASKS
%extract iceberg distributions for all DEMs without a manually-constructed dataset
disp('If any DEMs are new, adjust them for tides so sea level is at 0 m...');
close all; drawnow;
for p = 1:length(mats)

    %specify output location and file names
    cd([root_dir,'/',site_abbrev,'/DEMs/']);
    DEM_name = [site_abbrev,'-',DEMmat_dates(p,:),'_melange-DEM.mat']; 
    outputberg_name = [site_abbrev,'-',DEMmat_dates(p,:),'_melange-DEMfilled.mat'];
    
    %loop through DEMs & create a filled DEM with sea level at 0 m if one
    %does not already exist or if you specify (with filled_flag) that you
    %want to over-write the existing filled DEM
    if ismember(string(DEMmat_dates(p,:)),filledDEM_dates) == 0 || filled_flag == 1
        disp(['date ',num2str(p),' of ',num2str(length(mats)),' = ',DEMmat_dates(p,:)]);
        
        %load the existing filled DEM if specified
        if filled_flag == 1
            disp('loading gap-filled DEM');
            cd([root_dir,'/',site_abbrev,'/DEMs/']);
            load(outputberg_name); %load the gap-filled DEM
        end

        %adjust sea level to 0 m then fill holes in the DEM
        disp('loading raw DEM');
        cd([root_dir,'/',site_abbrev,'/DEMs/']);
        load(DEM_name); %load the masked orthometric melange elevations in Z structure
        disp('DEM loaded');
        
        %fill DEM gaps
        [Z,data_mask,gap_mask] = sl_correction(Z,WV_sigma); % perform sea level correction
        
        %save the DEM
        M.DEM.x = Z.x; M.DEM.y = Z.y; M.DEM.z = Z.z.adjusted;
        Z.z = rmfield(Z.z,'adjusted'); %only save necessary data to the elevation distribution data file
        save(DEM_name,'Z','-v7.3'); %SAVE
        
        %save fjord masks
        M.mask.fjord = Z.fjord.DEM_mask;
        M.mask.blunders = Z.melange.blunder_mask+gap_mask; M.mask.blunders(M.mask.blunders>1) = 1;
        M.mask.DEM = data_mask; clear data_mask;
        save(outputberg_name,'M','-v7.3');
        clear Z Y;
        
        %fill-in NaNs in the melange DEM
        disp('Filling NaNs in the melange DEM');
        fjord_elevs = M.mask.fjord.*M.DEM.z; fjord_elevs(fjord_elevs<0) = 0; fjord_elevs(M.mask.fjord==0) = 0;
        fjord_nans = sum(sum(isnan(fjord_elevs)));
        if fjord_nans < 100000
            infill_iters = 50;
        else %reduce inpaintn iterations to prevent computer from crashing (uses too much memory)
            infill_iters = 20;
        end
        fjord_filled = inpaintn(fjord_elevs,infill_iters); clear fjord_elevs;
        M.DEM.z_filled = fjord_filled; M.DEM.z_filled(M.mask.fjord~=1) = NaN;
        M.DEM.z_filled(M.mask.blunders==1) = NaN; M.DEM.z_filled = single(M.DEM.z_filled);
        M.DEM.z_filled(M.DEM.z_filled<0) = 0;
        M.DEM = rmfield(M.DEM,'z'); %M.DEM.z(M.DEM.z<0) = 0;
        save(outputberg_name,'M','-v7.3'); %SAVE
        
%         %if there is an "iceberg-data" matfile, which is leftover from an
%         %old processing pipeline, add it to the filled DEM then delete
%         for j = 1:size(iceberg_dates,1)
%             if contains(string(DEMmat_dates(p,:)),iceberg_dates(j,:))
%                 load(iceberg_mats(j).name);
%                 save(outputberg_name,'M','m','-v7.3');
%                 %move to an old file folder instead of deleting
%                 if exist('old-icebergs') == 0
%                     mkdir('old-icebergs')
%                 end
%                 movefile(iceberg_mats(j).name,['old-icebergs/',iceberg_mats(j).name]);
%                 
%                 %                 recycle('on'); %send to recycling, not permanently delete
%                 %                 delete_file = ['delete(''',iceberg_mats(j).name,''')']; eval(delete_file);
%             end
%         end
        
        %extract melange elevations
        melange = M.DEM.z_filled;
        melange(isnan(M.DEM.z_filled)) = 0;
        melange(melange<0) = 0; melange(melange>Hmax/(917/(1026-917))) = NaN;
        melange(M.mask.DEM==0) = NaN;

        %plot the melange DEM
        close(gcf);
        figure; h = histogram(melange(~isnan(melange)),min(melange(~isnan(melange))):1:max(melange(~isnan(melange)))); %use the elevation histogram to set colormap range
        cmax = h.BinEdges(find(cumsum(h.Values)>=0.98*sum(h.Values),1,'first')+1); elevcont_cmap = colormap(parula(max(threshold_range-adjuster)));
        clear h; close(gcf);
        figDEM = figure; set(figDEM,'position',[550 100 1000 500]);
        imagesc(M.DEM.x,M.DEM.y,melange); axis xy equal;
        melange_cmap = cmocean('thermal',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap);
        hold on; DEMax = gca;
        set(gca,'clim',[0 16]);  %set(gca,'clim',[0 cmax]);
        cbar = colorbar; cbar.Label.String  = 'elevation (m a.s.l.)';
        if find(nansum(melange)>0,1,'first')-50 > 1; mel_xlims = [find(nansum(melange)>0,1,'first')-50]; else; mel_xlims = 1; end
        if find(nansum(melange)>0,1,'last')+50 < length(M.DEM.x); mel_xlims = [mel_xlims find(nansum(melange)>0,1,'last')+50]; else; mel_xlims = [mel_xlims length(M.DEM.x)]; end
        if find(nansum(melange,2)>0,1,'first')-50 > 1; mel_ylims = [find(nansum(melange,2)>0,1,'first')-50]; else; mel_ylims = 1; end
        if find(nansum(melange,2)>0,1,'last')+50 < length(M.DEM.y); mel_ylims = [mel_ylims find(nansum(melange,2)>0,1,'last')+50]; else; mel_ylims = [mel_ylims length(M.DEM.y)]; end
        set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]); xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        title([DEMmat_dates(p,1:4),'/',DEMmat_dates(p,5:6),'/',DEMmat_dates(p,7:8)],'fontsize',16); grid on; drawnow;

        %uncomment next two lines if you want to delete the HUGE iceberg distribution files for the manual delineations
        %             recycle('on'); %send to recycling, not permanently delete
        %             delete([site_abbrev,'-',DEMmat_dates(p,:),'_melange-distribution.mat'])

        %clear a bunch of variables & close all figures
        clear Z *_mask G ib term* x_* y_* xi yi Zl* zjord* E...
            a ans bergH bergL bergW berg_W berg_length cmax fig*...
            fjord_nans fjord_filled in on infill_iters max_dim mel_*lims save_*...
            spurious_vals stats sub* z_adjusted z_mad z_median zfjord_filled...
            whalf berg_count termbuff_* m *grid *_sub zmin*;
        %             close all;

        %update list of filled DEM dates
        filled_DEMs = dir([site_abbrev,'*melange-DEMfilled.mat']);
        for i = 1:length(filled_DEMs)
            filledDEM_dates(i,:) = filled_DEMs(i).name(matfile_daterefs);
        end
        disp('Advancing...');

    end
    clear M;
end
disp("Done creating hole-less melange DEMs");

close all;

%% CHECK THAT MELANGE DEMS DO NOT CONTAIN SMALL REGIONS MISSED BY MANUAL MASKING
if mask_check == 1 
    disp('Manual check of melange masks...');
    %double-check the quality of the manual masks, removing small regions that were missed
    for p = 1:length(filled_DEMs)
        %load the data
        disp(['date ',num2str(p),' of ',num2str(length(filled_DEMs)),' = ',filledDEM_dates(p,:)]);
        DEM_name = [site_abbrev,'-',filledDEM_dates(p,:),'_melange-DEMfilled.mat'];
        load([root_dir,'/',site_abbrev,'/DEMs/',DEM_name]);
        [DEM_xgrid,DEM_ygrid] = meshgrid(M.DEM.x,M.DEM.y);

        %plot the map
        figDEM = figure; set(figDEM,'position',[50 100 500 1000]);
        sub1 = subplot(3,1,1);
        imagesc(M.DEM.x,M.DEM.y,M.DEM.z_filled); axis xy equal;
        melange_cmap = cmocean('thermal',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap);
        hold on; DEMax = gca;
        set(gca,'clim',[0 16]);  %set(gca,'clim',[0 cmax]);
        cbar = colorbar; cbar.Label.String  = 'elevation (m a.s.l.)';
        set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]); xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        title([filledDEM_dates(p,1:4),'/',filledDEM_dates(p,5:6),'/',filledDEM_dates(p,7:8)],'fontsize',16); grid on; drawnow;

        %make a dummy mask & plot
        BW_DEM = ones(size(M.DEM.z_filled));
        BW_DEM(isnan(M.DEM.z_filled)) = 0;
        BW_DEM(M.mask.DEM==0) = 0;
        figure(figDEM); sub2 = subplot(3,1,2);
        imagesc(M.DEM.x,M.DEM.y,M.mask.DEM); axis xy equal; colormap(gca,gray); hold on; DEMax = gca; grid on;
        set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]); xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        drawnow;

        %find any small unmasked bad or non-melange regions and mask them out
        stats = regionprops(logical(BW_DEM),'Area','PixelIdxList');
        for j = 1:size(stats,1); areas(j) = stats(j).Area; end
        [~,idx] = sort(areas,"descend");
        stats(idx(1)).PixelIdxList = NaN;
        for j = 1:size(stats,2)
            if ~isnan(stats(j).PixelIdxList)
                BW_DEM(stats(j).PixelIdxList) = 0;
            end
        end
        % pause
        M.mask.DEM = M.mask.DEM.*BW_DEM;

        %plot the mask
        figure(figDEM); subplot(sub2);
        imagesc(M.DEM.x,M.DEM.y,M.mask.DEM); axis xy equal; colormap(gca,gray); hold on; DEMax = gca;
        for j = 1:size(stats,2)
            if ~isnan(stats(j).PixelIdxList)
                plot(DEM_xgrid(stats(j).PixelIdxList),DEM_ygrid(stats(j).PixelIdxList),'.r'); hold on;
            end
        end
        set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]); xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);

        %manually mask-out any remaining splotches of data that were
        %somehow skipped (e.g., near fjord wall but connected to melange)
        % blunder_question = questdlg('Mask out any remaining "blunders"?',...
        %     'Blunder ID','1) Yes!','2) No!','1) Yes!');
        % switch blunder_question
        %     case '1) Yes!'
        %         figure(figDEM); subplot(sub1);
        %         disp('Click on UL & LR corners of a box bounding the anomalous elevations in the DEM subplot to zoom in'); % Upper left, lower right.
        %         [a] = ginput(2);
        %         set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
        %         drawnow;
        %         figure(figDEM); subplot(sub2);
        %         set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
        %         drawnow;
        %         disp('Draw a polygon on the mask to update it');
        %         [anom_zmask,xm,ym] = roipoly; anom_zmask = double(~anom_zmask);
        %         BW_DEM = anom_zmask.*BW_DEM;
        % 
        %         %update the plot
        %         imagesc(M.DEM.x,M.DEM.y,M.mask.DEM.*BW_DEM); axis xy equal; colormap(gca,gray); hold on; DEMax = gca;
        %         fill(xm,ym,'r'); hold on;
        %         set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]); xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        %         set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        %         xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        % 
        %         %update the mask
        %         M.mask.DEM = M.mask.DEM.*BW_DEM;
        %         clear anom_zmask xm ym;
        %     case '2) No!'
        %         disp('Good mask!')
        % end
        % clear blunder_question;

        %plot to check masking worked
        melange = M.DEM.z_filled;
        melange(isnan(M.DEM.z_filled)) = 0;
        melange(melange<0) = 0; melange(melange>Hmax/(917/(1026-917))) = NaN;
        melange(M.mask.DEM==0) = NaN;

        %update the figure
        figure(figDEM);
        sub3 = subplot(3,1,3);
        imagesc(M.DEM.x,M.DEM.y,melange); axis xy equal;
        melange_cmap = cmocean('thermal',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap);
        hold on; DEMax = gca; grid on;
        set(gca,'clim',[0 16]);  %set(gca,'clim',[0 cmax]);
        cbar = colorbar; cbar.Label.String  = 'elevation (m a.s.l.)';
        set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]); xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        drawnow;

        % %decide if you want to save the updated mask
        % resave_answer = questdlg('Save the updated mask?',...
        %     'mask update','1) Yes','2) No','2) No');
        % switch resave_answer
        %     case '1) Yes'
        %         save([root_dir,'/',site_abbrev,'/DEMs/',DEM_name],'M','-v7.3'); %SAVE
        %         disp('Updated the mask... moving on')
        %     case '2) No'
        %         disp('Did NOT update the mask... moving on')
        % end
        clear M DEM_name outputberg_name BW*DEM melange stats idx areas;
        % clear resave_answer;
        close(figDEM); drawnow;
    end
end

%% DELETE DEMS WITH HOLES TO REDUCE SPACE USE
% cd([root_dir,'/',site_abbrev,'/DEMs/']);
% for p = 1:length(mats)
%     recycle('on'); %send to recycling, not permanently delete
%     delete(mats(p).name);
% end

%% UPDATE FILE LIST
cd([root_dir,'/',site_abbrev,'/DEMs/']);
DEM_mats = dir([site_abbrev,'*_melange-DEMfilled.mat']); DEM_dates = ''; %DEMs
for p = 1:length(DEM_mats)
    DEM_dates(p,:) = DEM_mats(p).name(matfile_daterefs);
end
close all; drawnow;

%% EXTRACT ELEVATION TRANSECTS & AN ELEVATION PROFILE

%check if the elevation transects already exist (you've run the code
%before) and, if they do, prompt the user to decide if they want to re-run
%this section of code in order to update them
if exist([output_dir,site_abbrev,'_transects_elevations.csv']) == 2
    transect_answer = questdlg('Do you want to update the elevation transect CSV file?',...
        'transect update','1) Yes','2) No','2) No');
    switch transect_answer
        case '1) Yes'
            grab_profiles = 1; %if the transects have not been saved to a CSV file, extract them
        case '2) No'
            grab_profiles = 0; %skip elevation transect extraction
    end
else
    grab_profiles = 1; %if the transects have not been saved to a CSV file, extract them
end

%extract the profiles or move on
if grab_profiles == 1
    %initialize table creation
    y = []; x = []; column_names = ["Northing (m)","Easting (m)"];
    for j = 1:length(XF)
        hxf(j).x = XF(j).X'; hxf(j).y = XF(j).Y';
        y = [y; hxf(j).y; NaN];
        x = [x; hxf(j).x; NaN];
    end
    %check formatting of variables (must be column vectors)
    if size(AF.X,1) == 1; AF.X = AF.X'; AF.Y = AF.Y'; end
    %create tables to save elevation transects as CSVs
    Taf = table(AF.Y,AF.X); Taf.Properties.VariableNames = column_names;
    Txf = table(y,x); Txf.Properties.VariableNames = column_names;
    
    disp('Extracting along & across fjord elevation profiles...');
    for p = 1:length(DEM_mats)
        %load the data
        disp(['date ',num2str(p),' of ',num2str(length(DEM_mats)),' = ',DEM_dates(p,:)]);
        DEM_name = [site_abbrev,'-',DEM_dates(p,:),'_melange-DEMfilled.mat'];
        load([root_dir,'/',site_abbrev,'/DEMs/',DEM_name]);
        %     disp('Data loaded.');
        
        %fill data holes and masked regions with NaNs
        melange = M.DEM.z_filled;
        melange(melange<0) = NaN; melange(melange>Hmax/(rho_i/(rho_sw-rho_i))) = NaN;
        melange(M.mask.DEM==0) = NaN;
        
        %     %plot the DEM with transects if you want to check data coverage
        %     figure; set(gcf,'position',[50 50 800 800]);
        %     imagesc(M.DEM.x,M.DEM.y,melange); axis xy equal; hold on;
        %     melange_cmap = cmocean('thermal',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap);
        %     set(gca,'clim',[0 16],'fontsize',16); cbar = colorbar('fontsize',16); cbar.Label.String = 'elevation (m a.s.l.)';
        %     set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]);
        %     xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        %     set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        %     xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        %     plot(AF.X,AF.Y,'-k','linewidth',2);
        %     for j = 1:length(XF); plot(XF(j).X,XF(j).Y,'--k','linewidth',1.5); end
        
        %interpolate elevations to the centerline & add to the table
        [DEMx,DEMy] = meshgrid(M.DEM.x,M.DEM.y);
        h(1).z(:,p) = interp2(DEMx,DEMy,melange,AF.X,AF.Y);
        Taf.(['Elevation_',DEM_dates(p,:),' (m)']) = h(1).z(:,p);
        writetable(Taf,[output_dir,site_abbrev,'_centerline_elevations.csv']);
        
        %extract profiles and add to the data table: like a shapefile, each
        %transect is separated by a NaN for x,y,z values
        z = [];
        for j = 1:length(XF)
            h(j+1).z(:,p) = interp2(DEMx,DEMy,melange,XF(j).X',XF(j).Y');
            z = [z; h(j+1).z(:,p); NaN];
        end
        Txf.(['Elevation_',DEM_dates(p,:),' (m)']) = z;
        writetable(Txf,[output_dir,site_abbrev,'_transects_elevations.csv']);
        
        %clear variables and advance to the next date
        disp('Elevation CSVs updated & data added to "h" structure. Moving on to the next date.')
        clear M melange z;
        
    end
    clear Taf Txf;
    
    %load the Landsat reference image used to create the centerline & transects
    ims = dir([im_dir,'L*.TIF']);
    for j = 1:length(ims)
        if contains(ims(j).name,'B8')
            ref_image = [ims(j).folder,'/',ims(j).name];
        end
    end
    clear im_dir ims;
    %load the panchromatic Landsat scene & create a glacier/fjord mask
    [I,R] = readgeoraster(ref_image);
    %             info = geotiffinfo(ref_image);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I);
    clear I R;

    %crop the image to adjust brightnesses appropriately
    xlims = [find(im.x<=min(melmask.uncropped.x),1,'last'), find(im.x>=max(melmask.uncropped.x),1,'first')];
    ylims = [find(im.y>=max(melmask.uncropped.y),1,'last'), find(im.y<=min(melmask.uncropped.y),1,'first')];
    im_subset = im.z(min(ylims):max(ylims),min(xlims):max(xlims));
    im_subset = im_subset./max(max(im_subset));
    
    %plot the elevation time series as a quality check
    elev_fig = figure; set(gcf,'position',[50 50 2000 600]);
    sub1 = subplot(1,3,1); sub2 = subplot(1,3,2); sub3 = subplot(1,3,3);
    date_cmap = cmocean('solar',size(h(1).z,2)+1); tran_cmap = cmocean('deep',length(XF)+1);
    subplot(sub1);
    imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; hold on;
    plot(melmask.uncropped.x,melmask.uncropped.y,'-','color','w','linewidth',2); hold on;
    af_start = find(~isnan(nanmean(h(1).z,2)) == 1,1,'first');
    plot(AF.X(af_start:end),AF.Y(af_start:end),'-','color','k','linewidth',2); hold on;
    for j = 1:length(XF); plot(XF(j).X,XF(j).Y,'-','color',tran_cmap(j+1,:),'linewidth',2); hold on; end
    set(gca,'xlim',[min(AF.X(af_start:end))-2000,max(AF.X(af_start:end))+2000],...
        'ylim',[min(AF.Y(af_start:end))-2000,max(AF.Y(af_start:end))+2000],'fontsize',14);
    xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
    set(gca,'xticklabels',round(xticks/1000),'yticklabels',round(yticks/1000));
    xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16); drawnow;
    pos = get(sub1,'position'); set(sub1,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]);
    %plot the centerline profile time series
    subplot(sub2);
    af_dist(1) = 0;
    for k = af_start+1:length(AF.X)
        af_dist = [af_dist;af_dist(end)+sqrt((AF.X(k)-AF.X(k-1)).^2 + (AF.Y(k)-AF.Y(k-1)).^2)];
    end
    for j = 1:size(h(1).z,2); plot(af_dist,smooth(h(1).z(af_start:end,j),100),'-','color',date_cmap(j,:),'linewidth',2); hold on; end
    leg = legend(DEM_dates); leg.Location = 'west outside';
    if length(DEM_dates) > 30; leg.NumColumns = 2; end
    set(gca,'fontsize',14); grid on;
    xlabel('Distance along centerline (m)','fontsize',16);
    ylabel('Elevation (m a.s.l.)','fontsize',16);
    %plot the time-averaged transects
    subplot(sub3);
    for j = 1:length(XF)
        xf_dist(1) = 0;
        for k = 2:length(XF(j).X)
            xf_dist(k) = xf_dist(k-1)+sqrt((XF(j).X(k)-XF(j).X(k-1)).^2 + (XF(j).Y(k)-XF(j).Y(k-1)).^2);
        end
        plot(xf_dist,smooth(nanmedian(h(j+1).z,2),100),'-','color',tran_cmap(j+1,:),'linewidth',2); hold on;
        clear xf_dist;
    end
    set(gca,'fontsize',14); grid on;
    xlabel('Distance across fjord (m)','fontsize',16);
    ylabel('Elevation (m a.s.l.)','fontsize',16);
    pos = get(sub3,'position'); set(sub3,'position',[pos(1)+0.05 pos(2) pos(3) pos(4)]);
    set(sub2,'position',[0.45 pos(2) pos(3) pos(4)]);
    drawnow;
    saveas(elev_fig,[output_dir,site_abbrev,'-elevation-profiles.png'],'png');
    clear h;
    
    disp(['Completed elevation transect extraction! CSV at ',output_dir,site_abbrev,'_transects_elevations.csv']);
else
    disp(['Skipping elevation transect extraction b/c CSV already exists at ',output_dir,site_abbrev,'_transects_elevations.csv']);
end
clear grab_profiles;


%% check for overlapping transects and temporarily modify the melange mask to block transect overlap

%identify overlaps
for j = 1:length(XF)
    for k = 1:length(XF)
        if k ~= j
            [xi,~] = polyxpoly([XF(j).X(1),XF(j).X(end-1)],[XF(j).Y(1),XF(j).Y(end-1)],[XF(k).X(1),XF(k).X(end-1)],[XF(k).Y(1),XF(k).Y(end-1)]);
            if ~isempty(xi)
                transect_crosses(j,k) = k;
            else
                transect_crosses(j,k) = NaN;
            end
        else
            transect_crosses(j,k) = NaN;
        end
    end
end

%create a temporary melange mask that is modified as needed to provide a
%barrier for otherwise overlapping transects
if sum(~isnan(transect_crosses),'all') > 0
    transect_overlap = 1; %flag that the transects overlap and the melange mask is modifed
    
    %find the middle of neighboring intersections
    for j = 1:length(XF)
        if ismember(j-1,transect_crosses(j,:)) || ismember(j+1,transect_crosses(j,:))
            %         neighbor_crosses(j) = sum([ismember(j-1,transect_crosses(j,:)),ismember(j+1,transect_crosses(j,:))]);
            neighbor_crosses(j) = 1;
        else
            neighbor_crosses(j) = 0;
        end
    end
    mid_crosses = round(nanmean(find(neighbor_crosses==1)));
    
    %determine what side of the centerline the intersections occur on (assuming
    %it is just one side!) then find the segment of the fjord polygon that the
    %transect intersects on that side
    j = mid_crosses;
    for k = 1:length(XF)
        if k ~= j
            if ~isempty(polyxpoly([XF(j).X(1),XF(j).X(end-1)],[XF(j).Y(1),XF(j).Y(end-1)],[XF(k).X(1),XF(k).X(end-1)],[XF(k).Y(1),XF(k).Y(end-1)]))
                [xi(k),yi(k)] = polyxpoly([XF(j).X(1),XF(j).X(end-1)],[XF(j).Y(1),XF(j).Y(end-1)],[XF(k).X(1),XF(k).X(end-1)],[XF(k).Y(1),XF(k).Y(end-1)]);
            else
                xi(k) = NaN; yi(k) = NaN;
            end
        else
            xi(k) = NaN; yi(k) = NaN;
        end
    end
    %extend the transect so it intersects the melange outline
    dx1 = mode(diff(XF(j).X)); dy1 = mode(diff(XF(j).Y));
    if size(XF(j).X,1) == 1
        xd = [XF(j).X(1)-dx1,XF(j).X(1:end-1),XF(j).X(end-1)+dx1];
        yd = [XF(j).Y(1)-dy1,XF(j).Y(1:end-1),XF(j).Y(end-1)+dy1];
    else
        xd = [XF(j).X(1)-dx1;XF(j).X(1:end-1);XF(j).X(end-1)+dx1];
        yd = [XF(j).Y(1)-dy1;XF(j).Y(1:end-1);XF(j).Y(end-1)+dy1];
    end
    %find the intersections for the transect and the outline
    [xs,ys,is] = polyxpoly(xd,yd,melmask.uncropped.x,melmask.uncropped.y);
    %calculate the distance between each transect intersection point with the
    %middle of the intersecting transects and each melange edge
    dists(1,:) = sqrt((xs(1)-xi).^2 + (ys(1)-yi).^2);
    dists(2,:) = sqrt((xs(2)-xi).^2 + (ys(2)-yi).^2);
    %find which melange edge is closer to the transect intersections
    ind = find(dists == min(dists,[],'all')); [row,~] = ind2sub(size(dists),ind);
    
    %create a temp melange mask & add the part of the line that has
    %intersections to the shape
    melpoly_x = melmask.uncropped.x(1:is(row,2)); melpoly_y = melmask.uncropped.y(1:is(row,2));
    %only include the transect to the center-most intersection
    melpoly_x = [melpoly_x; xs(row); xi(find(dists(row,:) == max(dists(row,:)))); xs(row); melmask.uncropped.x(is(row,2)+1:end)];
    melpoly_y = [melpoly_y; ys(row); yi(find(dists(row,:) == max(dists(row,:)))); ys(row); melmask.uncropped.y(is(row,2)+1:end)];
    
else
    transect_overlap = 0;
    
    melpoly_x = melmask.uncropped.x;
    melpoly_y = melmask.uncropped.y;
end
clear *_crosses xi yi xd yd xs ys is dists ind row;

%% EXTRACT SIZE DISTRIBUTIONS

%convert elevation distributions into iceberg size distributions   
%note: this could be merged with the loop above to increase efficiency
disp('Extracting iceberg size distributions...');
for p = 1:length(DEM_mats)
    %load the data
    disp(['date ',num2str(p),' of ',num2str(length(DEM_mats)),' = ',DEM_dates(p,:)]);
    DEM_name = [site_abbrev,'-',DEM_dates(p,:),'_melange-DEMfilled.mat']; 
    load([root_dir,'/',site_abbrev,'/DEMs/',DEM_name]);
    disp('Data loaded.');
        
    %extract the distribution of elevations from the melange
    [M_xgrid,M_ygrid] = meshgrid(M.DEM.x,M.DEM.y);
    xvec = reshape(M_xgrid,1,[]); yvec = reshape(M_ygrid,1,[]); xy = [xvec; yvec];
    melange = M.DEM.z_filled;
    melange(isnan(M.DEM.z_filled)) = NaN;
    melange(melange<0) = 0; melange(melange>Hmax/(rho_i/(rho_sw-rho_i))) = NaN;
    melange(M.mask.DEM==0) = NaN;
%     figureA = figure; set(gcf,'position',[850 100 500 350]); title([site_abbrev,' ',DEM_dates(p,:),' iceberg size distribution']);
    [h.Values,h.BinEdges] = histcounts(melange(~isnan(melange)),0:1:ceil(Hmax/(rho_i/(rho_sw-rho_i)))); %hold on;
    bin_zno = double(h.Values); bin_zo = double((h.BinEdges(1:end-1) + h.BinEdges(2:end))/2);
    
    %convert to iceberg size distributions
    berg_W = (rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL.*bin_zo; Asurfs = (1/4)*pi*berg_W.^2; %calculate corresponding surface areas (m^2)
    Asurf_edges = (pi/4)*(rho_sw/(rho_sw-rho_i))^2*ARcomp.best.autoALL^2*h.BinEdges.^2; dA = Asurf_edges(2:end)-Asurf_edges(1:end-1); %calculate the bin widths in surface area units (m^2)
    berg_pixels = (pi*(0.5*(bin_zo.*((rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL))).^2)./(abs(M.DEM.x(1)-M.DEM.x(2)).^2);
    pixel_nos = bin_zno; berg_nos = pixel_nos./berg_pixels; pixels_per_bergclass = berg_pixels;
    m.melange.Asurfs = Asurfs; m.melange.bergs = berg_nos; m.melange.binwidth = dA; %add to iceberg size dataset
    %export the data as a table in a CSV
    T=table(berg_nos',Asurfs',dA');
    column_names = ["Count (unitless)", "Area (m^2)", "Area Binwidth (m^2)"];
    T.Properties.VariableNames = column_names; 
    writetable(T,[output_dir,site_abbrev,'-',num2str(DEM_dates(p,:)),'-iceberg-distribution.csv']);
    clear berg_* *_nos Asurf* dA T;
    
    %size distribution using the minimum aspect ratio: yields a smaller range of surface areas for icebergs
    berg_W = (rho_sw/(rho_sw-rho_i))*ARcomp.range.autoALL(1).*bin_zo; Asurfs = (1/4)*pi*berg_W.^2; 
    Asurf_edges = (pi/4)*(rho_sw/(rho_sw-rho_i))^2*ARcomp.range.autoALL(1)^2*h.BinEdges.^2; dA = Asurf_edges(2:end)-Asurf_edges(1:end-1);
    berg_pixels = (pi*(0.5*(bin_zo.*((rho_sw/(rho_sw-rho_i))*ARcomp.range.autoALL(1)))).^2)./(abs(M.DEM.x(1)-M.DEM.x(2)).^2);
    pixel_nos = bin_zno; berg_nos = pixel_nos./berg_pixels;
    m.melange.Asurfs_range(1,:) = Asurfs; m.melange.bergs_range(1,:) = berg_nos; m.melange.binwidth_range(1,:) = dA; %add to iceberg size dataset
    clear berg_* *_nos Asurf* dA;
    %size distribution using the maximum aspect ratio: yields a larger range of surface areas for icebergs
    berg_W = (rho_sw/(rho_sw-rho_i))*ARcomp.range.autoALL(2).*bin_zo; Asurfs = (1/4)*pi*berg_W.^2;
    Asurf_edges = (pi/4)*(rho_sw/(rho_sw-rho_i))^2*ARcomp.range.autoALL(2)^2*h.BinEdges.^2; dA = Asurf_edges(2:end)-Asurf_edges(1:end-1);
    berg_pixels = (pi*(0.5*(bin_zo.*((rho_sw/(rho_sw-rho_i))*ARcomp.range.autoALL(2)))).^2)./(abs(M.DEM.x(1)-M.DEM.x(2)).^2);
    pixel_nos = bin_zno; berg_nos = pixel_nos./berg_pixels;
    m.melange.Asurfs_range(2,:) = Asurfs; m.melange.bergs_range(2,:) = berg_nos; m.melange.binwidth_range(2,:) = dA; %add to iceberg size dataset
    clear berg_* *_nos Asurf* dA h; 
%     close(figureA); drawnow;
    save([root_dir,'/',site_abbrev,'/DEMs/',DEM_name],'M','m','-v7.3'); %iceberg size info
    
    %plot a 3-panel figure with (a) the melange DEM, 
    %(b) melange elevation histogram, and (c) the
    %full-melange iceberg size distribution in log-log space
    figure; set(gcf,'position',[50 50 500 1000]);
    sub1 = subplot(3,1,1); sub2 = subplot(3,1,2); sub3 = subplot(3,1,3);
    subplot(sub1);
    imagesc(M.DEM.x,M.DEM.y,melange); axis xy equal; hold on;
    melange_cmap = cmocean('thermal',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap);
    set(gca,'clim',[0 16],'fontsize',16); cbar = colorbar('fontsize',16); cbar.Label.String = 'elevation (m a.s.l.)';
    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]); 
    xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
    set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
    xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
    text(max(melmask.uncropped.x)-0.95*(max(melmask.uncropped.x)-min(melmask.uncropped.x)),min(melmask.uncropped.y)+0.925*(max(melmask.uncropped.y)-min(melmask.uncropped.y)),'a) melange DEM','fontsize',16); drawnow;
    drawnow;
    subplot(sub2);
    h = histogram(melange(~isnan(melange)),min(melange(~isnan(melange))):1:max(melange(~isnan(melange)))); hold on;
    xlabel('Melange elevation (m a.s.l.)','fontsize',16); ylabel('Pixel count','fontsize',16);
    set(gca,'fontsize',16,'YScale','log'); h.FaceColor = 'b'; h.FaceAlpha = 1; h.LineWidth = 1;
    xlims = get(sub2,'xlim'); ylims = get(sub2,'ylim');
    text(0.025*max(xlims),min(ylims)+0.28*(max(ylims)-min(ylims)),'b) melange elevation distribution','fontsize',16); drawnow;
    clear h; drawnow;
    subplot(sub3);
    fill([m.melange.Asurfs_range(1,m.melange.bergs_range(1,:)~=0) fliplr(m.melange.Asurfs_range(2,m.melange.bergs_range(2,:)~=0))],...
        [m.melange.bergs_range(1,m.melange.bergs_range(1,:)~=0) fliplr(m.melange.bergs_range(2,m.melange.bergs_range(2,:)~=0))],'b','EdgeColor','none','FaceAlpha',0.5); hold on;
    plot(m.melange.Asurfs(~isnan(m.melange.bergs)),m.melange.bergs(~isnan(m.melange.bergs)),'xk'); hold on;
%     leg = legend('observed distribution','fontsize',18);
    set(gca,'fontsize',16,'XScale','log','YScale','log',...
        'ytick',[10^-6,10^-4,0.01,1,100,10000,1000000],'yticklabel',[10^-6,10^-4,0.01,1,100,10000,1000000],...
        'xtick',[100,1000,10000,100000,1000000],'xticklabel',[100,1000,10000,100000,1000000]); %grid on; 
    xlims = get(sub3,'xlim'); ylims = get(sub3,'ylim');
    xlabel('Iceberg surface area (m^2)','fontsize',16); ylabel('Iceberg count','fontsize',16);
    text(min(xlims)+0.00003*(max(xlims)-min(xlims)),min(ylims)+0.28*(max(ylims)-min(ylims)),'c) iceberg size distribution','fontsize',16); drawnow; clear xlims ylims;
    drawnow;
    saveas(gcf,[output_dir,site_abbrev,'-',num2str(DEM_dates(p,:)),'_melange-distribution_subplots.png'],'png');
    disp('Saved subplots & textfile of size distributions ');
    
    %compile info for all dates
    ib(p).Asurf = m.melange.Asurfs; ib(p).bergs = m.melange.bergs; 
    
    %repeat size distribution extraction but for a subset of the melange:
    %iteratively move from the seaward limit to the terminus
    Tsub = table(m.melange.Asurfs',m.melange.binwidth');
    column_names = ["Area (m^2)", "Area Binwidth (m^2)"];
    Tsub.Properties.VariableNames = column_names; 
    for j = 1:length(XF)-1
        disp(['Subsetting melange: ',num2str(j),' of ',num2str(length(XF)-1)]);
        %note: slightly extend the transects to make sure that Matlab finds
        %an intersection on both end with the melange outline 
        
        %find the centroid coordinates for the polygon
        [xij,yij,ij] = polyxpoly([XF(j).X,XF(j).X(end-1)],[XF(j).Y,XF(j).Y(end-1)],AF.X,AF.Y);
        [xijp,yijp,ijp] = polyxpoly([XF(j+1).X,XF(j+1).X(end-1)],[XF(j+1).Y,XF(j+1).Y(end-1)],AF.X,AF.Y);
        centroid_coords = [xij yij; xijp yijp];
        
        %first transect: extend & identify intersections
        dx1 = mode(diff(XF(j).X)); dy1 = mode(diff(XF(j).Y)); 
        if size(XF(j).X,1) == 1
            XF(j).X = [XF(j).X(1)-dx1,XF(j).X(1:end-1),XF(j).X(end-1)+dx1];
            XF(j).Y = [XF(j).Y(1)-dy1,XF(j).Y(1:end-1),XF(j).Y(end-1)+dy1];
        else
            XF(j).X = [XF(j).X(1)-dx1;XF(j).X(1:end-1);XF(j).X(end-1)+dx1];
            XF(j).Y = [XF(j).Y(1)-dy1;XF(j).Y(1:end-1);XF(j).Y(end-1)+dy1];
        end
        [xi1,yi1,i1] = polyxpoly(XF(j).X,XF(j).Y,melpoly_x,melpoly_y);
        %second transect: extend & identify intersections
        dx2 = mode(diff(XF(j+1).X)); dy2 = mode(diff(XF(j+1).Y)); 
        if size(XF(j+1).X,1) == 1
            XF(j+1).X = [XF(j+1).X(1)-dx2,XF(j+1).X(1:end-1),XF(j+1).X(end-1)+dx2];
            XF(j+1).Y = [XF(j+1).Y(1)-dy1,XF(j+1).Y(1:end-1),XF(j+1).Y(end-1)+dy2];
        else
            XF(j+1).X = [XF(j+1).X(1)-dx2;XF(j+1).X(1:end-1);XF(j+1).X(end-1)+dx2];
            XF(j+1).Y = [XF(j+1).Y(1)-dy2;XF(j+1).Y(1:end-1);XF(j+1).Y(end-1)+dy2];
        end
        [xi2,yi2,i2] = polyxpoly(XF(j+1).X,XF(j+1).Y,melpoly_x,melpoly_y);
        
        %if there are >2 intersections because the fjord curved, filter out
        %the ones farthest from the centerline... the 'first' & 'last' are
        %somewhat arbitrary but prevent double identification of non-unique
        %intersections where a hinge in the polygons doubles-back
        if size(i1,1) > 2 %first transect
            di1 = i1(:,1)-ij(:,1);
            ref1a = find(di1==max(di1(di1<0)),1,'last'); ref1b = find(di1==min(di1(di1>0)),1,'first'); 
            xi1_temp = xi1([ref1a,ref1b]); yi1_temp = yi1([ref1a,ref1b]); i1_temp = i1([ref1a,ref1b],:); 
            clear xi1 yi1 i1;
            xi1 = xi1_temp; yi1 = yi1_temp; i1 = i1_temp; 
            clear di1 ref1* *i1_temp;
        end
        if size(i2,1) > 2 %second transect
            di2 = i2(:,1)-ijp(:,1);
            ref2a = find(di2==max(di2(di2<0)),1,'last'); ref2b = find(di2==min(di2(di2>0)),1,'first'); 
            xi2_temp = xi2([ref2a,ref2b]); yi2_temp = yi2([ref2a,ref2b]); i2_temp = i2([ref2a,ref2b],:); 
            clear xi2 yi2 i2;
            xi2 = xi2_temp; yi2 = yi2_temp; i2 = i2_temp; 
            clear di2 ref2* *i2_temp;
        end
        
        %note: the i1,i2 values should give the indices of the fjord outline segment that is intersected by each transect
        [~,nntei] = min(sqrt(sum(([xi1(1),yi1(1)] - [xi2,yi2]).^2,2))); %find the nearest end vertex for the 2nd transect 
        if i2(nntei,2) > i1(1,2)
            mel1_inds = [i1(1,2)+1:1:i2(nntei,2)]; %fjord mask indices increase between transect ends
            %confirm that the fjord sidewall length is reasonable with a
            %positive increment between indices (a huge distance means that
            %the indices span the starting point of the melange mask)
            side_coords = [xi1(1) yi1(1); melpoly_x(mel1_inds) melpoly_y(mel1_inds); xi2(nntei) yi2(nntei)];
            if sum(sqrt(sum(diff(side_coords).^2,2))) > 2*sum(sqrt(sum(diff([xi1(1) yi1(1); xi2(nntei) yi2(nntei)]).^2,2))) %apply an arbitrary threshold based on centroid distances for transects
                clear mel1_inds;
                mel1_inds = [i1(1,2):-1:1, length(melpoly_x):-1:i2(nntei,2)+1]; %fjord mask indices wrap around melange mask start/end
            end
        elseif i2(nntei,2) < i1(1,2)
            mel1_inds = [i1(1,2):-1:i2(nntei,2)+1];
            side_coords = [xi1(1) yi1(1); melpoly_x(mel1_inds) melpoly_y(mel1_inds); xi2(nntei) yi2(nntei)];
            if sum(sqrt(sum(diff(side_coords).^2,2))) > 2*sum(sqrt(sum(diff([xi1(1) yi1(1); xi2(nntei) yi2(nntei)]).^2,2))) %apply an arbitrary threshold based on centroid distances for transects
                mel1_inds = [i1(1,2)+1:1:length(melpoly_x), 1:1:i2(nntei,2)]; %fjord mask indices wrap around melange mask start/end
            end
        elseif i2(nntei,2) == i1(1,2)
            mel1_inds = [];
        else
            error('no melange mask coordinates identified!');
        end
        clear side_coords;
        %repeat but for the other transect ends
        [~,nntef] = min(sqrt(sum(([xi1(2),yi1(2)] - [xi2,yi2]).^2,2))); %find the nearest end vertex for the 2nd transect 
        if i2(nntef,2) < i1(2,2)
            mel2_inds = [i2(nntef,2)+1:1:i1(2,2)]; %fjord mask indices increase between transect ends
            %confirm that the fjord sidewall length is reasonable with a
            %positive increment between indices (a huge distance means that
            %the indices span the starting point of the melange mask)
            side_coords = [xi2(nntef) yi2(nntef); melpoly_x(mel2_inds) melpoly_y(mel2_inds); xi1(2) yi1(2)];
            if sum(sqrt(sum(diff(side_coords).^2,2))) > 2*sum(sqrt(sum(diff([xi1(2) yi1(2); xi2(nntef) yi2(nntef)]).^2,2))) %apply an arbitrary threshold based on centroid distances for transects
                clear mel2_inds;
                mel2_inds = [i2(nntef,2):-1:1, length(melpoly_x):-1:i1(2,2)+1]; %fjord mask indices wrap around melange mask start/end
            end
        elseif i2(nntef,2) > i1(2,2)
            mel2_inds = [i2(nntef,2):-1:i1(2,2)+1];
            side_coords = [xi2(nntef) yi2(nntef); melpoly_x(mel2_inds) melpoly_y(mel2_inds); xi1(2) yi1(2)];
            if sum(sqrt(sum(diff(side_coords).^2,2))) > 2*sum(sqrt(sum(diff([xi1(2) yi1(2); xi2(nntef) yi2(nntef)]).^2,2))) %apply an arbitrary threshold based on centroid distances for transects
                mel2_inds = [i2(nntef,2)+1:1:length(melpoly_x), 1:1:i1(2,2)]; %fjord mask indices wrap around melange mask start/end
            end
        elseif i2(nntef,2) == i1(2,2)
            mel2_inds = [];
        else
            error('no melange mask coordinates identified!');
        end
        clear side_coords;
        %put the full polygon together to subset the melange
        if size(XF(j).X,1) ~= size(melpoly_x,1) && size(XF(j).X,1)==1
            melsubset_polyx = [xi1(1), melpoly_x(mel1_inds)', xi2(nntei), xi2(nntef), melpoly_x(mel2_inds)', xi1(2), xi1(1)];
            melsubset_polyy = [yi1(1), melpoly_y(mel1_inds)', yi2(nntei), yi2(nntef), melpoly_y(mel2_inds)', yi1(2), yi1(1)];
        elseif size(XF(j).X,1) ~= size(melpoly_x,1) && size(melpoly_x,1)==1
            melsubset_polyx = [xi1(1), melpoly_x(mel1_inds), xi2(nntei), xi2(nntef), melpoly_x(mel2_inds), xi1(2), xi1(1)];
            melsubset_polyy = [yi1(1), melpoly_y(mel1_inds), yi2(nntei), yi2(nntef), melpoly_y(mel2_inds), yi1(2), yi1(1)];
        else
            error('Unexpected polyline shapes');
        end

        %if the melange mask protruded to the side for real (say where a
        %tributary glacier enters the system) and the index adjustment to
        %account for wrapping around the ends forced the mask to go the
        %other way, the centroid should not be in the mask & you need to
        %force re-do it (Ex: Midgard Glacier, SE Greenland = MGG)
        [stat] = inpoly2(nanmean(centroid_coords),[melsubset_polyx',melsubset_polyy']);
        if stat == 0
            clear mel*_inds melsubset_polyx melsubset_polyy;
            %redo first set of side coordinates
            if i2(nntei,2) > i1(1,2)
                mel1_inds = [i1(1,2)+1:1:i2(nntei,2)];
            elseif i2(nntei,2) < i1(1,2)
                mel1_inds = [i1(1,2):-1:i2(nntei,2)+1];
            elseif i2(nntei,2) == i1(1,2)
                mel1_inds = [];
            else
                error('no melange mask coordinates identified!');
            end
            clear side_coords;
            %redo second set of side coordinates
            if i2(nntef,2) < i1(2,2)
                mel2_inds = [i2(nntef,2)+1:1:i1(2,2)]; %fjord mask indices increase between transect ends
            elseif i2(nntef,2) > i1(2,2)
                mel2_inds = [i2(nntef,2):-1:i1(2,2)+1];
            elseif i2(nntef,2) == i1(2,2)
                mel2_inds = [];
            else
                error('no melange mask coordinates identified!');
            end
            clear side_coords;
            %put the full polygon together to subset the melange
            if size(XF(j).X,1) ~= size(melpoly_x,1) && size(XF(j).X,1)==1
                melsubset_polyx = [xi1(1), melpoly_x(mel1_inds)', xi2(nntei), xi2(nntef), melpoly_x(mel2_inds)', xi1(2), xi1(1)];
                melsubset_polyy = [yi1(1), melpoly_y(mel1_inds)', yi2(nntei), yi2(nntef), melpoly_y(mel2_inds)', yi1(2), yi1(1)];
            elseif size(XF(j).X,1) ~= size(melpoly_x,1) && size(melpoly_x,1)==1
                melsubset_polyx = [xi1(1), melpoly_x(mel1_inds), xi2(nntei), xi2(nntef), melpoly_x(mel2_inds), xi1(2), xi1(1)];
                melsubset_polyy = [yi1(1), melpoly_y(mel1_inds), yi2(nntei), yi2(nntef), melpoly_y(mel2_inds), yi1(2), yi1(1)];
            else
                error('Unexpected polyline shapes');
            end
        end
        clear stat;
        clear xi* yi* i1 i2 mel*_inds;
        
        %plot the melange mask subset as a quality check
        if p == 1 %only plot for the first DEM to check that it is masking correctly
            if j == 1
                subfig = figure; set(gcf,'position',[50 50 500 500]);
                imagesc(M.DEM.x,M.DEM.y,melange); axis xy equal; hold on;
                melange_cmap = cmocean('thermal',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap);
                set(gca,'clim',[0 16],'fontsize',16); cbar = colorbar('fontsize',16); cbar.Label.String = 'elevation (m a.s.l.)';
                set(gca,'xlim',[min(melmask.dated(p).x) max(melmask.dated(p).x)],'ylim',[min(melmask.dated(p).y) max(melmask.dated(p).y)]);
                xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
                set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
                xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
                plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',2); hold on;
            end
            figure(subfig);
            plot(XF(j).X,XF(j).Y,'-b','linewidth',2); hold on;
            plot(XF(j+1).X,XF(j+1).Y,'-b','linewidth',2); hold on;
            plot(melsubset_polyx,melsubset_polyy,'--c','linewidth',2); hold on;
            drawnow;
        end
            
        %find the melange elevations inside the polygon
%             disp(['identifying DEM pixels in subset...']);
%             tic
        [stat] = inpoly2(xy',[melsubset_polyx',melsubset_polyy']);
        in = reshape(stat,size(melange));
%             toc
        [hsub.Values,hsub.BinEdges] = histcounts(melange(in),0:1:ceil(Hmax/(rho_i/(rho_sw-rho_i)))); hold on;
        binsub_zno = double(hsub.Values); %binsub_zo = double((hsub.BinEdges(1:end-1) + hsub.BinEdges(2:end))/2);
        bergsub_nos = binsub_zno./pixels_per_bergclass;
        %uncomment the next 7 lines to create test maps for
        %identification of DEM pixels inside the subsetted melange polygon
        % figure; imagesc(x_subset,y_subset,in); axis xy equal; hold on;
        % mask_cmap = colormap('gray'); mask_cmap = flipud(mask_cmap); colormap(gca,mask_cmap); colorbar;
        % set(gca,'xlim',[min(melmask.dated(p).x) max(melmask.dated(p).x)],'ylim',[min(melmask.dated(p).y) max(melmask.dated(p).y)]);
        % xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        % set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        % xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        % plot(melmask.uncropped.x,melmask.uncropped.y,'-b','linewidth',2); hold on;
        % plot(melsubset_polyx,melsubset_polyy,'--c','linewidth',2); hold on;

        %add to the table if the melange area is big enough
        pixel_area = abs((M.DEM.x(2)-M.DEM.x(1)).*(M.DEM.y(2)-M.DEM.y(1)));
        centroid = mean(centroid_coords);
        disp(['percent data coverage = ',num2str(100*sum(~isnan(melange(in)))./sum(in(~isnan(in))))]);
        disp(['area of data coverage = ',num2str(sum(~isnan(melange(in))).*pixel_area),' m^2']);
        % if sum(~isnan(melange(in)))./sum(in(~isnan(in))) > 0.25 %if more than 25% of the subset has data, save it
        if sum(~isnan(melange(in))).*pixel_area > (transect_spacer*1000) %if the subset area with data covers an average of >1 km-wide, save it
            Tsub.([num2str(round(centroid(1))),'E,',num2str(round(centroid(2))),'N']) = bergsub_nos';
            writetable(Tsub,[output_dir,site_abbrev,'-',num2str(DEM_dates(p,:)),'-iceberg-distribution-subsets.csv']);
            disp('subset data saved to CSV');
        else
            Tsub.([num2str(round(centroid(1))),'E,',num2str(round(centroid(2))),'N']) = NaN(size(hsub.Values))';
            writetable(Tsub,[output_dir,site_abbrev,'-',num2str(DEM_dates(p,:)),'-iceberg-distribution-subsets.csv']);
            disp('not enough elevation data in the subset... writing NaNs to CSV');
        end
        
        clear melsubset* in hsub binsub* bergsub* pixelsub* centroid* *_inds *_subset stat;
    end
    
    %close the figure as needed and clear date-specific variables
    if p == 1; close(subfig); end
    clear bin* melange M m T Tsub;
end
close all; drawnow;

%plot size distributions color-coded according to date
DEMyrs = str2num(DEM_dates(:,1:4)); DEMmos = str2num(DEM_dates(:,5:6));
sizedist_fig = figure; set(sizedist_fig,'position',[1050 50 700 200*ceil((max(DEMyrs)-min(DEMyrs)+1)/2)+100]);
month_cmap = cmocean('phase',12);
% [unique_yrs,unique_yrrefs] = unique(DEMyrs); 
%create subplots
for j = 1:(max(DEMyrs)-min(DEMyrs)+1)
    subplot(ceil((max(DEMyrs)-min(DEMyrs)+1)/2),2,j);
end
%plot data
for p = 1:numel(DEM_mats)
    ib(p).mel_area = sum(ib(p).Asurf.*ib(p).bergs);
    subplot(ceil((max(DEMyrs)-min(DEMyrs)+1)/2),2,DEMyrs(p)-min(DEMyrs)+1);
%         disp(['melange area = ',num2str(ib(p).mel_area./10^6),'km^2']);
%     if ib(p).mel_area >= 175e6
        loglog(ib(p).Asurf(ib(p).bergs>0.01),ib(p).bergs(ib(p).bergs>0.01),'-','color',month_cmap(DEMmos(p),:),'linewidth',1.5); hold on;
%     end
    
end
%label subplots
for j = 1:(max(DEMyrs)-min(DEMyrs)+1)
    subplot(ceil((max(DEMyrs)-min(DEMyrs)+1)/2),2,j); axpos = get(gca,'position');
    set(gca,'position',[axpos(1)-0.05*(1-mod(j,2))-0.02 axpos(2)+0.07*(ceil(j/2)./ceil((max(DEMyrs)-min(DEMyrs)+1)/2)) axpos(3) axpos(4)]); axpos = get(gca,'position'); 
    loglog(ib(1).Asurf(ib(1).bergs>0.01),ib(1).bergs(ib(1).bergs>0.01),'-','color','none','linewidth',1.5); hold on;
    set(gca,'ytick',[0.01,1,100,10000,1000000],'yticklabel',[],...
        'xtick',[100,1000,10000,100000,1000000],'xticklabel',[]); set(gca,'fontsize',16);
    text(70,0.05,[alphabet(j),') ',num2str(min(DEMyrs)+j-1)],'fontsize',16);
    %label axes
    if j == (2*ceil((max(DEMyrs)-min(DEMyrs)+1)/2))-1
        set(gca,'yticklabel',[0.01,1,100,10000,1000000],'xticklabel',[100,1000,10000,100000,1000000],...
            'xticklabelrotation',45); 
        xlabel('Iceberg surface area (m^2)','fontsize',16); ylabel('Iceberg count','fontsize',16);
    end
    %add colorbar
    if j == max(DEMyrs)-min(DEMyrs)+1
        for k = 1:12; pl(k) = loglog(ib(end).Asurf(1),ib(end).bergs(1),'-','color',month_cmap(k,:),'linewidth',1.5); hold on; end
        leg=legend(pl,{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'},'location','eastoutside',...
            'orientation','horizontal'); legpos = get(leg,'position');
        set(leg,'fontsize',16,'NumColumns',1); set(leg,'position',[0.885 0.55-0.5*legpos(4) 0.1 legpos(4)]);
    end
    %enlarge subplots
    set(gca,'position',[axpos(1) axpos(2) 1.1*axpos(3) 1.1*axpos(4)]);
end
saveas(sizedist_fig,[output_dir,site_abbrev,'-iceberg-size-distribution-timeseries-subplots.png'],'png');
disp(['Finished extracting iceberg size distributions from DEMs for ',site_abbrev,'!']);
disp(['... full-melange size distributions are saved as ',output_dir,site_abbrev,'-yyyymmdd-iceberg-distribution.csv']);
disp(['... subsetted size distributions are saved as ',output_dir,site_abbrev,'-yyyymmdd-iceberg-distribution-subsets.csv']);


%% compile all the time-stamped size distributions into a single CSV
berg_dists = dir([output_dir,site_abbrev,'*-iceberg-distribution.csv']);

T = table;
for p = 1:length(berg_dists)
    %read the CSV
    Ttemp = readtable([output_dir,berg_dists(p).name]);
    Tarray = table2array(Ttemp);
    msg = warning('query','last'); wid = msg.identifier; warning('off',wid);

    %add data to the site-compiled table
    if p == 1
        T.('Area (m^2)') = Tarray(:,2); T.('Area Binwidth (m^2)') = Tarray(:,3);
    end
    T.(['Count_',berg_dists(p).name(5:12),' (unitless)']) = Tarray(:,1);
    
    clear Ttemp Tarray;
end
clear berg_dists;

%export the data compiled for the site as a single table
writetable(T,[output_dir,site_abbrev,'-iceberg-distribution-timeseries.csv']);
clear T;

end

