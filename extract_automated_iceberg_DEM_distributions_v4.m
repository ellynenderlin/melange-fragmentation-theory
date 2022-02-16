function extract_automated_iceberg_DEM_distributions_v4(root_path,glacier_abbrev,output_path)
%------------------------------------
% Description: Construct iceberg size distributions from
% automatically-extracted elevation distributions
%
%may need to add some paths (depending on permissions) before running: 
%run in the command window when specifying root_path, glacier_abbrev,output_path
% addpath('/users/ellynenderlin/mfiles','/users/ellynenderlin/mfiles/general',...
%     '/users/ellynenderlin/mfiles/Greenland-calving/icebergs',...
%     '/users/ellynenderlin/mfiles/Greenland-calving/icebergs/automated-iceberg-mapping',...
%     '/users/ellynenderlin/mfiles/Greenland-calving/icebergs/fragtheory_data_fit');
%
%------------------------------------------
% %% set manually
% clear all; 
% 
% basepath='/users/ellynenderlin/iceberg-fragmentation/';
% root_path = basepath;
% glacier_abbrev = 'KL';
% output_path = basepath;

%% specify variables
rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)
WV_sigma = 2.9; %DEM uncertainty in meters (used when filtering bad data)
z_inc = ceil(WV_sigma); %increment for melange elevation bins (min ~= WV_sigma; controls number of iceberg size classes)
threshold_range = [z_inc:1:12]; %test range for elevation thresholds to delineate icebergs using elevation contours from automated_manual_iceberg_aspects_comparison.m (needed for plot colors)
adjuster = z_inc-1; %scalar used to assign outputs of thresholding to a structure starting at index 1
Wmax = 1000; %maximum likely iceberg width (based on qualitative inspection of images during terminus mapping)
Hmax = 800; %threshold thickness in meters that you do not expect icebergs to exceed (grounding line thickness is a good proxy)
ARcomp.best.autoALL = 2; % iceberg aspect ratio 
ARcomp.range.autoALL = [1.7, 2.3]; % range in iceberg aspect ratio

%load data
cd_to_output = ['cd ',output_path,'/',glacier_abbrev,'/']; eval(cd_to_output);
load_mask = ['load ',glacier_abbrev,'-melange-masks.mat']; eval(load_mask); %created using create_melange_masks_v2
mask_dates = '';
for i = 1:length(melmask.dated)
    mask_dates(i,:) = melmask.dated(i).datestring;
end

%extract decimal dates for DEMs
mats = dir([glacier_abbrev,'*melange-DEM.mat']);
leap_doys = [31 29 31 30 31 30 31 31 30 31 30 31]; norm_doys = [31 28 31 30 31 30 31 31 30 31 30 31];
leap_cumdoys = cumsum(leap_doys); leap_cumdoys = [0 leap_cumdoys(1:11)];
norm_cumdoys = cumsum(norm_doys); norm_cumdoys = [0 norm_cumdoys(1:11)];
for i = 1:length(mats)
    DEMmat_dates(i,:) = mats(i).name(4:11);
    if mod(str2num(DEMmat_dates(1:4)),4) ~= 0
        DEM_decidate(i,:) = str2num(DEMmat_dates(i,1:4)) + (norm_cumdoys(str2num(DEMmat_dates(i,5:6)))+str2num(DEMmat_dates(i,7:8)))./sum(norm_doys);
    else
        DEM_decidate(i,:) = str2num(DEMmat_dates(i,1:4)) + (leap_cumdoys(str2num(DEMmat_dates(i,5:6)))+str2num(DEMmat_dates(i,7:8)))./sum(leap_doys);
    end
end
filled_DEMs = dir([glacier_abbrev,'*melange-DEMfilled.mat']);
for i = 1:length(filled_DEMs)
    filledDEM_dates(i,:) = filled_DEMs(i).name(4:11);
end

%identify existing melange iceberg distribution datasets
cd_to_glacier = ['cd ''',root_path,'/',glacier_abbrev,'''']; eval(cd_to_glacier);
melange_mats = dir([glacier_abbrev,'*_melange-distribution.mat']); melange_dates = ''; %manual iceberg sizes
for i = 1:length(melange_mats)
    melange_dates(i,:) = melange_mats(i).name(4:11);
end
iceberg_mats = dir([glacier_abbrev,'*_iceberg-data.mat']); iceberg_dates = ''; %automated iceberg sizes
for i = 1:length(iceberg_mats)
    iceberg_dates(i,:) = iceberg_mats(i).name(4:11);
end

alphabet = 'abcdefghijklmnopqrstuvwxyz';
%% PROCESS MELANGE MASKS
%extract iceberg distributions for all DEMs without a manually-constructed dataset
disp(['Extracting iceberg distributions for ',glacier_abbrev,'...']);
close all; drawnow;
for p = 1:length(mats)

    %specify output location and file names
    cd_to_output = ['cd ',output_path,'/',glacier_abbrev,'/']; eval(cd_to_output);
    DEM_name = [glacier_abbrev,'-',DEMmat_dates(p,:),'_melange-DEM.mat']; save_DEM = ['save(''',DEM_name,''',''Z'',''-v7.3'')']; %raw & intermediate elevation data
    %         filledDEM_name = [glacier_abbrev,'-',DEMmat_dates(p,:),'_melange-DEMfilled.mat']; %filled, orthorectified DEM
    %         outputberg_name = [glacier_abbrev,'-',DEMmat_dates(p,:),'_iceberg-data.mat'];
    outputberg_name = [glacier_abbrev,'-',DEMmat_dates(p,:),'_melange-DEMfilled.mat'];
    %identify the original & filled DEMs
    for k = 1:size(filledDEM_dates,1)
        filledflag(k) = contains(string(DEMmat_dates(p,:)),filledDEM_dates(k,:));
    end
    
    %loop through DEMs WITHOUT manual aspect ratios & extract elevation
    %distributions (already did this step for the dates with manual data)
    if isempty(strmatch(string(DEMmat_dates(p,:)),string(melange_dates))) == 1
        disp(['date ',num2str(p),' of ',num2str(length(mats)),' = ',DEMmat_dates(p,:)]);
        %load the gap-filled DEM if already generated or the raw DEM if new
        if sum(filledflag)>0
            disp('loading gap-filled DEM');
            cd_to_output = ['cd ',output_path,'/',glacier_abbrev,'/']; eval(cd_to_output);
            load_DEM = ['load ',filled_DEMs(filledflag==1).name]; eval(load_DEM); %load the gap-filled DEM
            clear filledflag;
        end
        %check that filled DEM contains all the necessary data
        if exist('M') ~=1
            disp('loading raw DEM');
            cd_to_glacier = ['cd ''',root_path,'/',glacier_abbrev,'''']; eval(cd_to_glacier);
            load_DEM = ['load ',DEM_name]; eval(load_DEM); %load the masked orthometric melange elevations in Z structure
            disp('DEM loaded');
            
            %fill DEM gaps
            [Z,data_mask,gap_mask] = sl_correction(Z,WV_sigma); % perform sea level correction
            eval(save_DEM); %SAVE
            M.DEM.x = Z.x; M.DEM.y = Z.y; M.DEM.z = Z.z.adjusted; Z.z = rmfield(Z.z,'adjusted'); %only save necessary data to the elevation distribution data file
            eval(save_DEM); %SAVE
            
            %save fjord masks
            M.mask.fjord = Z.fjord.DEM_mask;
            M.mask.blunders = Z.melange.blunder_mask+gap_mask; M.mask.blunders(M.mask.blunders>1) = 1;
            M.mask.DEM = data_mask; clear data_mask;
            save_filledDEM = ['save(''',outputberg_name,''',''M'',''-v7.3'')']; eval(save_filledDEM);
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
            eval(save_filledDEM); %SAVE
        end
        
        %if there is an "iceberg-data" matfile, which is leftover from an
        %old processing pipeline, add it to the filled DEM then delete
        for j = 1:length(iceberg_dates)
            if contains(string(DEMmat_dates(p,:)),iceberg_dates(j,:))
                load_old_file = ['load ',iceberg_mats(j).name]; eval(load_old_file);
                save_filledDEM = ['save(''',outputberg_name,''',''M'',''m'',''-v7.3'')']; eval(save_filledDEM);
                %move to an old file folder instead of deleting
                if exist('old-icebergs') == 0
                    mkdir('old-icebergs')
                end
                movefile(iceberg_mats(j).name,['old-icebergs/',iceberg_mats(j).name]);
                
                %                 recycle('on'); %send to recycling, not permanently delete
                %                 delete_file = ['delete(''',iceberg_mats(j).name,''')']; eval(delete_file);
            end
        end
%         if exist('m','var')
%             save_icebergs = ['save(''',outputberg_name,''',''bergs'',''m'',''-v7.3'')']; %iceberg size info
%             disp('saving');
%         else
%             save_icebergs = ['save(''',outputberg_name,''',''bergs'',''-v7.3'')']; %iceberg size info
%             disp('saving');
%         end
        
        %check if icebergs have already been automatically delineated
        %TEMPORARILY REMOVED THIS CHECK WHILE TESTING
%         if isempty(strmatch(string(DEMmat_dates(p,:)),string(filledDEM_dates))) == 1            

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
            
%             %save figures
%             saveas(figDEM,[glacier_abbrev,'-',DEMmat_dates(p,:),'_melange-DEM.png'],'png'); 
            
            %uncomment next two lines if you want to delete the HUGE iceberg distribution files for the manual delineations
%             recycle('on'); %send to recycling, not permanently delete
%             delete([glacier_abbrev,'-',DEMmat_dates(p,:),'_melange-distribution.mat'])
            
          %clear a bunch of variables & close all figures
            clear Z *_mask G M ib term* x_* y_* xi yi Zl* zjord* E...
                a ans bergH bergL bergW berg_W berg_length cmax fig*...
                fjord_nans fjord_filled in on infill_iters max_dim mel_*lims save_*...
                spurious_vals stats sub* z_adjusted z_mad z_median zfjord_filled...
                whalf berg_count termbuff_* m *grid *_sub zmin*;
            close all;
            disp('Advancing');
%         end
        
    else
        disp('Iceberg data already extracted');
    end
end
disp("Done.");
%% UPDATE FILES AND PLOT ICEBERG SIZE DISTRIBUTIONS
cd_to_output = ['cd ',output_path,'/',glacier_abbrev,'/']; eval(cd_to_output);
DEM_mats = dir([glacier_abbrev,'*_melange-DEMfilled.mat']); DEM_dates = ''; %DEMs
for i = 1:length(DEM_mats)
    DEM_dates(i,:) = DEM_mats(i).name(4:11);
end
close all; drawnow;

%convert elevation distributions into iceberg size distributions   
%note: this could be merged with the loop above to increase efficiency
disp('Creating a time series of iceberg size distributions...');
for p = 1:length(DEM_mats)
    cd_to_output = ['cd ',output_path,'/',glacier_abbrev,'/']; eval(cd_to_output);
    
    %load the data
    disp(['date ',num2str(p),' of ',num2str(length(DEM_mats)),' = ',DEM_dates(p,:)]);
    DEM_name = [glacier_abbrev,'-',DEM_dates(p,:),'_melange-DEMfilled.mat']; load_DEM = ['load ',DEM_name]; eval(load_DEM);
    outputberg_name = [glacier_abbrev,'-',DEM_dates(p,:),'_melange-DEMfilled.mat'];
%     bergs_name = [glacier_abbrev,'-',DEM_dates(p,:),'_iceberg-data.mat']; load_bergs = ['load ',bergs_name]; eval(load_bergs);
%     outputberg_name = [glacier_abbrev,'-',DEM_dates(p,:),'_iceberg-data.mat']; save_icebergs = ['save(''',outputberg_name,''',''bergs'',''m'',''-v7.3'')']; %iceberg size info
    disp('Data loaded.');
        
    %extract the distribution of elevations from the melange
    melange = M.DEM.z_filled;
    melange(isnan(M.DEM.z_filled)) = 0;
    melange(melange<0) = 0; melange(melange>Hmax/(rho_i/(rho_sw-rho_i))) = NaN;
    melange(M.mask.DEM==0) = NaN;
    figureA = figure; set(gcf,'position',[850 100 500 350]); title([glacier_abbrev,' ',DEM_dates(p,:),' iceberg size distribution']);
    h = histogram(melange(~isnan(melange)),0:1:ceil(Hmax/(rho_i/(rho_sw-rho_i)))); hold on;
    bin_zno = double(h.Values); bin_zo = double((h.BinEdges(1:end-1) + h.BinEdges(2:end))/2);
    
    %convert to iceberg size distributions
    berg_W = (rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL.*bin_zo; Asurfs = (1/4)*pi*berg_W.^2; %calculate corresponding surface areas (m^2)
    Asurf_edges = (pi/4)*(rho_sw/(rho_sw-rho_i))^2*ARcomp.best.autoALL^2*h.BinEdges.^2; dA = Asurf_edges(2:end)-Asurf_edges(1:end-1); %calculate the bin widths in surface area units (m^2)
    berg_pixels = (pi*(0.5*(bin_zo.*((rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL))).^2)./(abs(M.DEM.x(1)-M.DEM.x(2)).^2);
    pixel_nos = bin_zno; berg_nos = pixel_nos./berg_pixels;
%     max_bin = find(berg_nos>0.1,1,'last'); %identify fractions of bergs <0.1berg
%     pixel_nos(max_bin+1) = sum(pixel_nos(max_bin+1:end)); pixel_nos(max_bin+2:end) = NaN; %merge fractions of bergs <0.1berg
    berg_nos = pixel_nos./berg_pixels;
    dlmwrite([glacier_abbrev,'-',num2str(DEM_dates(p,:)),'-iceberg-distribution.txt'],[berg_nos' Asurfs' dA'],'delimiter','\t','precision','%.1f'); %export the full melange size distribution data
    m.melange.Asurfs = Asurfs; m.melange.bergs = berg_nos; m.melange.binwidth = dA; %add to iceberg size dataset
    clear berg_* *_nos Asurf* dA;
    %size distribution using the minimum aspect ratio: yields a smaller range of surface areas for icebergs
    berg_W = (rho_sw/(rho_sw-rho_i))*ARcomp.range.autoALL(1).*bin_zo; Asurfs = (1/4)*pi*berg_W.^2; 
    Asurf_edges = (pi/4)*(rho_sw/(rho_sw-rho_i))^2*ARcomp.range.autoALL(1)^2*h.BinEdges.^2; dA = Asurf_edges(2:end)-Asurf_edges(1:end-1);
    berg_pixels = (pi*(0.5*(bin_zo.*((rho_sw/(rho_sw-rho_i))*ARcomp.range.autoALL(1)))).^2)./(abs(M.DEM.x(1)-M.DEM.x(2)).^2);
    pixel_nos = bin_zno; berg_nos = pixel_nos./berg_pixels;
    berg_nos = pixel_nos./berg_pixels;
    m.melange.Asurfs_range(1,:) = Asurfs; m.melange.bergs_range(1,:) = berg_nos; m.melange.binwidth_range(1,:) = dA; %add to iceberg size dataset
    clear berg_* *_nos Asurf* dA;
    %size distribution using the maximum aspect ratio: yields a larger range of surface areas for icebergs
    berg_W = (rho_sw/(rho_sw-rho_i))*ARcomp.range.autoALL(2).*bin_zo; Asurfs = (1/4)*pi*berg_W.^2;
    Asurf_edges = (pi/4)*(rho_sw/(rho_sw-rho_i))^2*ARcomp.range.autoALL(2)^2*h.BinEdges.^2; dA = Asurf_edges(2:end)-Asurf_edges(1:end-1);
    berg_pixels = (pi*(0.5*(bin_zo.*((rho_sw/(rho_sw-rho_i))*ARcomp.range.autoALL(2)))).^2)./(abs(M.DEM.x(1)-M.DEM.x(2)).^2);
    pixel_nos = bin_zno; berg_nos = pixel_nos./berg_pixels;
    berg_nos = pixel_nos./berg_pixels;
    m.melange.Asurfs_range(2,:) = Asurfs; m.melange.bergs_range(2,:) = berg_nos; m.melange.binwidth_range(2,:) = dA; %add to iceberg size dataset
    clear berg_* *_nos Asurf* dA h; close(figureA); drawnow;
%     bergs = m.melange.bergs; % ADDED BY JUKES - CONFIRM THIS IS OKAY  
%     outputberg_name = [glacier_abbrev,'-',DEMmat_dates(p,:),'_melange-DEMfilled.mat'];        
%     if exist('m','var')
%         save_icebergs = ['save(''',outputberg_name,''',''bergs'',''m'',''-v7.3'')']; %iceberg size info
%         disp('saving');
%     else
%         save_icebergs = ['save(''',outputberg_name,''',''bergs'',''-v7.3'')']; %iceberg size info
%         disp('saving');
%     end
    save_icebergs = ['save(''',outputberg_name,''',''M'',''m'',''-v7.3'')']; %iceberg size info
    eval(save_icebergs); %SAVE
    
    %plot a 3-panel figure with (a) the melange DEM, 
    %(b) melange elevation histogram, and (c) the
    %full-melange iceberg size distribution in log-log space
    figure; set(gcf,'position',[1000 500 500 1000]);
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
    save_fig = ['saveas(gcf,''',glacier_abbrev,'-',num2str(DEM_dates(p,:)),'_melange-distribution_subplots.png'',''png'');']; eval(save_fig);
    disp('Saved subplots & textfile of size distributions ');
    
    %compile info for all dates
    ib(p).Asurf = m.melange.Asurfs; ib(p).bergs = m.melange.bergs; 

    clear bin* melange M m;
end
    
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
saveas(sizedist_fig,[glacier_abbrev,'_iceberg-size-distribution-timeseries-subplots.png'],'png');
end

