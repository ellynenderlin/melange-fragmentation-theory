%%% Create plots & maps of GrIS melange data

clearvars; close all; warning off;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/inpoly2/');
addpath('/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-fragmentation-code/');

%specify paths
% root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-melt/';
root_dir = '/Volumes/Jokulhaup_5T/Greenland-melange/';
years = 2011:1:2023; yr_cmap = cmocean('matter',length(years)+1); yr_cmap = yr_cmap(2:end,:);
mo_cmap = cmocean('phase',12); close all;

%identify the site folders
cd(root_dir);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.') && length(sites(i).name) == 3
        sitenames = [sitenames; sites(i).name];
    end
end
%if you have a preferred order for the sites in plots, specify it here
geo_order = [3,14,11,4,10,8,2,13,12,6,9,1,7,5,15]; [refs,sort_ind] = sort(geo_order);
% sort_ind = [1:1:length(sitenames)]; %if you have no predefined site order, progress alphabetically

%option to reload the existing data if any data have been saved
if exist([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']) == 2
    reload = questdlg('Reload the existing centerline data & start from there?',...
        'data reload','1) Yes: reload','2) No: start fresh','1) Yes: reload');
    switch reload
        case '1) Yes: reload'
            load([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']);
            if length(MP) < size(sitenames,1)
                site_start = length(MP)+1;
                disp(['...restarting on site #',num2str(site_start),' (',sitenames(site_start,:),')']);
            else
                disp('data reloaded but dataset is fully processed');
                replot = questdlg('Have you commented-out data processing code to focus only on plotting?',...
                    'replot','1) Yes: replot','2) No!','1) Yes: replot');
                switch replot
                    case '1) Yes: replot'
                        site_start = 1;
                    case '2) No!'
                        disp('comment-out data processing code if you only want to replot, then rerun');
                        return
                end
            end
        case '2) No: start fresh'
            site_start = 1;
            MP = struct;
    end
end

%% loop through the folders & extract info
disp('Creating plots of elevation, velocity, and terminus position...');
for j = site_start:length(sitenames)
    disp(sitenames(j,:)); output_dir = [root_dir,sitenames(j,:),'/'];
    MP(j).name = sitenames(j,:);
    
    %navigate to the study site directory
    cd([root_dir,sitenames(j,:)]);
    LCdir = dir([root_dir,sitenames(j,:),'/LC*']); im_dir = [LCdir(1).folder,'/',LCdir(1).name,'/']; %Landsat 8 or 9 unzipped image directory for mapping
    
    %load the time-stamped melange masks and extract the approximate
    %terminus position as the intersection between the centerline and the
    %DEM-based melange outline
    load([MP(j).name,'-melange-masks.mat']); %load the melange mask file
    DEM_num = size(melmask.dated,2); term_trace = [];
    
    %load the shapefile of transect-centerline intersections used to
    %extract the velocity timeseries & the full centerline to precisely
    %pinpoint relative movement of the terminus position
    C = readtable([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'_centerline_2000m-interval.csv']);
    MP(j).V.X = C.Easting_m_; MP(j).V.Y = C.Northing_m_; clear C;
    clear C;
    C = shaperead([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'_centerline.shp']);
    centerline_dist = [0, cumsum(sqrt((C.X(2:end)-C.X(1:end-1)).^2 + (C.Y(2:end)-C.Y(1:end-1)).^2))]';
    for l = 1:length(MP(j).V.X)
        dists = sqrt((MP(j).V.X(l)-C.X).^2 + (MP(j).V.Y(l)-C.Y).^2);
        tran_dist(l) = centerline_dist(find(dists == min(dists)));
        clear dists;
    end
    
%     %identify which terminus traces are actually the edge of the DEM
%     %(straight lines) then make an overview map that excludes them
%     ims = dir([im_dir,'L*.TIF']);
%     for k = 1:length(ims)
%         if contains(ims(k).name,'B8')
%             ref_image = [ims(k).folder,'/',ims(k).name];
%         end
%     end
%     clear im_dir ims;
%     %load the panchromatic Landsat scene
%     [I,R] = readgeoraster(ref_image);
%     im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
%     im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
%     im.z = double(I);
%     clear I R;
%     %crop the image to adjust brightnesses appropriately
%     xlims = [find(im.x<=min(melmask.uncropped.x),1,'last'), find(im.x>=max(melmask.uncropped.x),1,'first')];
%     ylims = [find(im.y>=max(melmask.uncropped.y),1,'last'), find(im.y<=min(melmask.uncropped.y),1,'first')];
%     im_subset = im.z(min(ylims):max(ylims),min(xlims):max(xlims));
%     im_subset = im_subset./max(max(im_subset));
    
    %if you are just replotting the map, uncomment the lines immediately
    %below and then skip to the map plotting section of code
        for p = 1:DEM_num
            term_trace(p) = MP(j).T.qualflag(p);
            zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
            datest = char(MP(j).Z.date(p));
            yrs(p) = str2num(datest(1:4)); mos(p) = str2num(datest(5:6));
            clear datest;
        end
        term_ref = find(abs(zdate-2020.66) == min(abs(zdate(term_trace==1)-2020.66))); %use terminus delineation closest to Aug. 2020 as the centerline reference
        tran_reldist = MP(j).T.centerline(1,term_ref) - tran_dist;
        centerline_reldist = MP(j).T.centerline(1,term_ref) - centerline_dist;
    
%     %plot the image
%     temp_fig = figure; set(temp_fig,'position',[50 50 1200 1200]);
%     for p = 1:DEM_num
%         imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; drawnow; hold on;
%         set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)],'fontsize',18);
%         plot(C.X,C.Y,'.k'); hold on;
%         plot(melmask.dated(p).x,melmask.dated(p).y,'-m','linewidth',2); drawnow;
%         answer = questdlg('Is the terminus delineation good?',...
%             'terminus delineation','1) Yes: wiggly & good','2) No: DEM edge (straight)','3) No: wonky polygon','1) Yes: wiggly & good');
%         switch answer
%             case '1) Yes: wiggly & good'
%                 term_trace = [term_trace; 1]; melmask.dated(length(term_trace)).terminus = 1;
%                 MP(j).T.qualflag(p) = 1;
%             case '2) No: DEM edge (straight)'
%                 term_trace = [term_trace; 0]; melmask.dated(length(term_trace)).terminus = 1;
%                 MP(j).T.qualflag(p) = 0;
%             case '3) No: wonky polygon'
%                 removed_flag = fix_individual_melange_masks(root_dir,sitenames(j,:),melmask,p);
%                 if strmatch(removed_flag,'removed')
%                     %removed the DEM from melmask so reload it
%                     load([MP(j).name,'-melange-masks.mat']);
%                     DEM_num = size(melmask.dated,2); p = p-1; %reset counters to account for removed data
%                 else
%                     term_trace = [term_trace; 1]; melmask.dated(length(term_trace)).terminus = 1;
%                     MP(j).T.qualflag(p) = 1;
%                 end
%         end
%         clear answer; cla;
%     end
%     close(temp_fig);
%     save([root_dir,sitenames(j,:),'/',sitenames(j,:),'-melange-masks.mat'],'melmask','-v7.3');
%     
    %open the CSV of time-stamped elevation transects
    T = readtable([sitenames(j,:),'_transects_elevations.csv']);
    T_headers = T.Properties.VariableNames;
    datestart = strfind(string(T_headers(3)),'20');
    for p = 3:size(T,2)
        MP(j).Z.date(1,p-2) = string(T_headers{p}(datestart:datestart+7));
    end
%     
%     %intersect the centerline with each time-stamped melange outline
%     for p = 1:size(melmask.dated,2)
%         zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
%         [xis,yis,iis] = polyxpoly(melmask.dated(p).x,melmask.dated(p).y,MP(j).V.X,MP(j).V.Y);
%         MP(j).T.X(1,p) = xis(end); MP(j).T.Y(1,p) = yis(end);
%         if term_trace(p) == 1 %terminus was mapped from the DEM
%             if length(xis)>1
%                 MP(j).T.centerline(1,p) = tran_dist(iis(end,2))+sqrt((MP(j).T.X(1,p)-MP(j).V.X(iis(end,2))).^2 + (MP(j).T.Y(1,p)-MP(j).V.Y(iis(end,2))).^2);
%             else
%                 MP(j).T.centerline(1,p) = tran_dist(iis(1,2))+sqrt((MP(j).T.X(1,p)-MP(j).V.X(iis(1,2))).^2 + (MP(j).T.Y(1,p)-MP(j).V.Y(iis(1,2))).^2);
%             end
%         else %terminus was cut-off in the DEM so don't record the centerline intersection
%             MP(j).T.centerline(1,p) = NaN;
%         end
%         clear xis yis iis;
%     end
%     term_ref = find(abs(zdate-2020.66) == min(abs(zdate(term_trace==1)-2020.66))); %use terminus delineation closest to Aug. 2020 as the centerline reference
%     tran_reldist = MP(j).T.centerline(1,term_ref) - tran_dist;
%     centerline_reldist = MP(j).T.centerline(1,term_ref) - centerline_dist;
%     
%     %create an overview map of all the dated melange masks
%     map_fig = figure; set(map_fig,'position',[850 50 800 600]);
%     imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; drawnow; hold on;
%     set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]);
%     %plot dummy lines to create a legend
%     first_term = find(term_trace == 1,1,'first');
%     for l = 1:size(yr_cmap,1)
%         pt(l) = plot(melmask.dated(first_term).x,melmask.dated(first_term).y,...
%             '-','color',yr_cmap(l,:),'linewidth',1.5); hold on;
%         leg_labels(l) = cellstr(num2str(years(l)));
%     end
%     %plot the real terminus traces
%     for p = 1:size(melmask.dated,2)
%         if term_trace(p) == 1
%             plot(melmask.dated(p).x,melmask.dated(p).y,'-','color',...
%                 yr_cmap(str2num(melmask.dated(p).datestring(1:4))-min(years)+1,:),'linewidth',1.5); hold on;
%             drawnow;
%         end
%     end
%     %plot the centerline
%     plot(C.X,C.Y,'.b'); hold on;
%     %replot the terminus trace used as the reference position for distance 
%     pt(size(yr_cmap,1)+1) = plot(melmask.dated(term_ref).x,melmask.dated(term_ref).y,'-g','linewidth',1.5); hold on;
%     leg_labels(size(yr_cmap,1)+1) = {'reference'};
%     %plot the uncropped melange mask
%     plot(melmask.uncropped.x,melmask.uncropped.y,'-b','linewidth',2); hold on;
%     %add transect locations
%     plot(MP(j).V.X,MP(j).V.Y,'.k','markersize',16,'linewidth',2); hold on;
%     xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
%     set(gca,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
%     xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
%     pos = get(gca,'position'); 
%     if range(xlims) > 1.25*range(ylims) %short and fat map so plot the legend below
%         map_leg = legend(pt,leg_labels,'Location','southoutside',...
%             'Orientation','horizontal','NumColumns',7);
%         leg_pos = get(map_leg,'position');
%         set(map_leg,'position',[leg_pos(1)+(0.5-mean([leg_pos(1) leg_pos(1)+leg_pos(3)])) 0.075 leg_pos(3) leg_pos(4)]);
%         if range(xlims)/range(ylims) < 1.5
%             set(gca,'position',[pos(1) 0.225 0.9*pos(3) 0.9*pos(4)]); drawnow;
%         else
%             set(gca,'position',[pos(1) 0.125 pos(3) pos(4)]); drawnow;
%         end
%     else %tall and thin map so plot the legend on the side
%         map_leg = legend(pt,leg_labels,'Location','eastoutside',...
%             'Orientation','vertical','NumColumns',1);
%         leg_pos = get(map_leg,'position');
%         if leg_pos(1)+0.05 > pos(1)-0.05+pos(3)
%             set(map_leg,'position',[leg_pos(1)+0.05 leg_pos(2) leg_pos(3) leg_pos(4)]);
%         else
%             set(map_leg,'position',[pos(1)-0.05+pos(3) leg_pos(2) leg_pos(3) leg_pos(4)]);
%         end
%         set(gca,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]); drawnow;
%     end
%     saveas(map_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'_site-map.png'],'png'); %save the image
%     
    %find the NaNs in the coordinate pairs to identify each transect
    coords = table2array(T(:,1:2));
    nan_inds = find(isnan(coords(:,1))==1);
    
    %extract the data
    tran_ind = 1; k = 1; zprofs = [];
    while k < size(T,1)
        if k == 1
            %extract full coordinates for transect
            MP(j).Z.Xrange(tran_ind,1:2) = coords([1,nan_inds(k)-1],2);
            MP(j).Z.Yrange(tran_ind,1:2) = coords([1,nan_inds(k)-1],1);

            %extract the centroid coordinates
            MP(j).Z.X(tran_ind,1) = nanmean(table2array(T(1:nan_inds(k)-1,2)));
            MP(j).Z.Y(tran_ind,1) = nanmean(table2array(T(1:nan_inds(k)-1,1)));

            %calculate the mean elevation for each date
            ztemp = table2array(T(1:nan_inds(k)-1,3:end)); ztemp(ztemp<3) = NaN;
            MP(j).Z.Zavg(tran_ind,:) = nanmean(ztemp); clear ztemp;

            %add to a temp matrix for plotting
            zprofs(tran_ind,:) = MP(j).Z.Zavg(tran_ind,:);

            tran_ind = tran_ind+1;
        elseif ismember(k,nan_inds)
            if ~ismember(k+1,nan_inds) %check for back-to-back NaNs
                %extract full coordinates for transect
                MP(j).Z.Xrange(tran_ind,1:2) = coords([k+1,nan_inds(find(nan_inds==k)+1)-1],2);
                MP(j).Z.Yrange(tran_ind,1:2) = coords([k+1,nan_inds(find(nan_inds==k)+1)-1],1);

                %extract the centroid coordinates
                MP(j).Z.X(tran_ind,1) = nanmean(table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,2)));
                MP(j).Z.Y(tran_ind,1) = nanmean(table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,1)));

                %calculate the median elevation for each date
                ztemp = table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,3:end)); ztemp(ztemp<3) = NaN;
                MP(j).Z.Zavg(tran_ind,:) = nanmean(ztemp); clear ztemp;

                %add to a temp matrix for plotting
                zprofs(tran_ind,:) = MP(j).Z.Zavg(tran_ind,:);

                tran_ind = tran_ind+1;
            end
        end
        k = k+1;
    end
    clear T T_headers datestart coords nan_inds;
% 
%     %load the velocity timeseries for the transect-centerline intersection
%     %points and plot a velocity profile with the closest mid-date to each
%     %elevation profile (if after 2013, when Landsat 8 was launched)
%     vel_pts = dir([root_dir,sitenames(j,:),'/velocities/']);
%     for i = 1:length(vel_pts)
%         if contains(vel_pts(i).name,'velocity')
%             pt_ref = str2num(vel_pts(i).name(end-5:end-4));
% 
%             %read the file
%             V = readtable([root_dir,sitenames(j,:),'/velocities/',vel_pts(i).name]);
% 
%             %filter out all the velocities based on temporal resolution
%             short_dts = find(V.days_dt<60); %get rid of all velocities with coarse temporal resolution
%             vel_dates = V.mid_date(short_dts); vel_dts = V.days_dt(short_dts);
%             vels = V.velocity_m_yr_(short_dts); vels(vels == 0) = NaN;
% 
%             %convert datetime to a decimal date
%             for k = 1:length(vel_dates)
%                 decidate(k) = convert_to_decimaldate(vel_dates(k),'yyyy-MM-dd HH:mm:ss.SSS');
%             end
% 
%             %create an average velocity profile using all velocity
%             %observations with a mid-date within 15 days of the DEM
%             for p = 1:length(MP(j).Z.date)
%                 if zdate(p) > 2013
%                     datediff = abs(zdate(p) - decidate);
%                     MP(j).V.date(pt_ref,p) = string(nanmean(datetime(vel_dates(find(abs(datediff) <= (15/365))),'Format','yyyyMMdd')));
%                     MP(j).V.dt(pt_ref,p) = nanmean(vel_dts(find(abs(datediff) <= (15/365))));
%                     MP(j).V.V(pt_ref,p) = nanmean(vels(find(abs(datediff) <= (15/365))));
%                     clear datediff;
%                 else
%                     MP(j).V.date(pt_ref,p) = NaN; MP(j).V.dt(pt_ref,p) = NaN;
%                     MP(j).V.V(pt_ref,p) = NaN;
%                 end
%             end
% 
%             clear pt_ref V short_dts vel_dates vel_dts vels decidate;
%         end
%     end

    %plot time-series of the width-averaged elevation profiles
    figure; set(gcf,'position',[50 50 1200 700]);
    subZ_yr = subplot(2,3,1); subZ_mo = subplot(2,3,4); 
    mean_prof = nanmean(zprofs,2); 
    seaward_idx = find(~isnan(mean_prof)==1,1,'first');
    inland_idx = find(~isnan(mean_prof)==1,1,'last');
    for p = 1:size(zprofs,2)
        %add dummy lines for the legend
        if p == 1
            for l = 1:size(yr_cmap,1)
                subplot(subZ_yr);
                py(l) = plot(tran_reldist(seaward_idx:end)',zprofs(seaward_idx:end,p),'-','color',yr_cmap(l,:),'linewidth',2); hold on;
            end
            for l = 1:size(mo_cmap,1)
                subplot(subZ_mo);
                pm(l) = plot(tran_reldist(seaward_idx:end)',zprofs(seaward_idx:end,p),'-','color',mo_cmap(l,:),'linewidth',2); hold on;
            end
        end
        
        %plot annual averages
        subplot(subZ_yr);
        MP(j).Z.Zavg(MP(j).Z.Zavg==0) = NaN; %remove mysterious zeros (need to look into their source!)
        for k = 1:length(years)
            yr_idx = find(yrs == years(k));
            plot(tran_reldist(seaward_idx:end)',nanmean(MP(j).Z.Zavg(seaward_idx:end,yr_idx),2),'-','color',yr_cmap(years(k)-min(years)+1,:),'linewidth',2); hold on; 
            clear yr_idx
        end
        
        %plot all data color-coded by season
        yr = str2num(MP(j).Z.date{p}(1:4)); mo = str2num(MP(j).Z.date{p}(5:6));
%         subplot(subZ_yr);
%         plot(tran_reldist(seaward_idx:end)',zprofs(seaward_idx:end,p),'-','color',yr_cmap(yr-min(years)+1,:),'linewidth',2); hold on; 
        subplot(subZ_mo);
        plot(tran_reldist(seaward_idx:end)',zprofs(seaward_idx:end,p),'-','color',mo_cmap(mo,:),'linewidth',2); hold on;
        clear yr mo;
    end
    subplot(subZ_yr);
    title([sitenames(j,:),' elevation profiles'])
    plot(tran_reldist(seaward_idx:end)',mean_prof(seaward_idx:end),'-','color','k','linewidth',3); hold on; 
    pos = get(subZ_yr,'position'); 
    leg_yr = legend(py,num2str(years')); leg_yr.Location = 'eastoutside';
    set(subZ_yr,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]); %shift plot back to pre-legend location
    set(gca,'fontsize',18); grid on; drawnow;
    subplot(subZ_mo); 
    plot(tran_reldist(seaward_idx:end)',mean_prof(seaward_idx:end),'-','color','k','linewidth',3); hold on; 
    pos = get(subZ_mo,'position'); 
    leg_mo = legend(pm,num2str([1:1:12]')); leg_mo.Location = 'eastoutside';
    set(subZ_mo,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]); %shift plot back to pre-legend location
    set(gca,'fontsize',18); grid on; drawnow;
    Zxlims = get(subZ_mo,'xlim');
    
    %create velocity profile plots
    subV_yr = subplot(2,3,2); subV_mo = subplot(2,3,5);
    for p = 1:length(MP(j).Z.date)
        %plot annual averages
        subplot(subV_yr);
        for k = 1:length(years)
            yr_idx = find(yrs == years(k));
            plot(tran_reldist(seaward_idx:end)',nanmean(MP(j).V.V(seaward_idx:end,yr_idx),2)./365,'-','color',yr_cmap(years(k)-min(years)+1,:),'linewidth',2); hold on; 
            clear yr_idx
        end
        
        %plot all data color-coded by season
        yr = str2num(MP(j).Z.date{p}(1:4)); mo = str2num(MP(j).Z.date{p}(5:6));
%         subplot(subV_yr);
%         plot(tran_reldist(seaward_idx:end)',MP(j).V.V(seaward_idx:end,p)./365,'-','color',yr_cmap(yr-min(years)+1,:),'linewidth',2); hold on;
        subplot(subV_mo);
        plot(tran_reldist(seaward_idx:end)',MP(j).V.V(seaward_idx:end,p)./365,'-','color',mo_cmap(mo,:),'linewidth',2); hold on;
        clear yr mo;
    end
    subplot(subV_yr); 
    plot(tran_reldist(seaward_idx:end)',nanmean(MP(j).V.V(seaward_idx:end,:),2)./365,'-','color','k','linewidth',3); hold on; 
    title([sitenames(j,:),' velocity profiles'])
    set(gca,'fontsize',18); grid on; drawnow;
    subplot(subV_mo); 
    plot(tran_reldist(seaward_idx:end)',nanmean(MP(j).V.V(seaward_idx:end,:),2)./365,'-','color','k','linewidth',3); hold on; 
    set(gca,'fontsize',18); grid on; drawnow;
    Vxlims = get(subV_mo,'xlim');
    
    %create a terminus position timeseries
%     term_fig = figure;
    subT = subplot(2,3,[3,6]);
    for p = 1:length(MP(j).Z.date)
%         if term_trace(p) == 1 %don't plot terminus delineations that are the DEM edge, not the true terminus
            mo = str2num(MP(j).Z.date{p}(5:6));
            plot(MP(j).T.centerline(1,term_ref)-MP(j).T.centerline(p),zdate(p),'x','color',mo_cmap(mo,:),'linewidth',2); hold on;
            clear mo;
%         end
    end
    set(subT,'xlim',[min(MP(j).T.centerline(1,term_ref)-MP(j).T.centerline),max(ceil(tran_reldist(seaward_idx:inland_idx)/1000)*1000)],'ylim',[min(years) max(years)]); 
    xlims = get(subT,'xlim'); xticks = get(subT,'xtick');
    set(gca,'xtick',xticks,'xticklabels',xticks/1000,'fontsize',18);
    xlabel('Distance from terminus (km)','fontsize',18); ylabel('Year','fontsize',18); 
    grid on; drawnow;
%     Txlims = get(gca,'xlim');
    
    %standardize the x-limits on the plots
%     set(subZ_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subZ_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
%     set(subV_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subV_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
    subplot(subZ_yr);
    set(subZ_yr,'xlim',xlims,'xtick',xticks,'xticklabels',[]); 
    pos = get(subZ_yr,'position'); set(subZ_yr,'position',[pos(1) pos(2)-0.05 pos(3) pos(4)+0.05]); %slightly stretch y-axis
    ylabel('Elevation (m)','fontsize',18); 
    subplot(subZ_mo);
    set(subZ_mo,'xlim',xlims,'xtick',xticks,'xticklabels',xticks/1000);
    pos = get(subZ_mo,'position'); set(subZ_mo,'position',[pos(1) pos(2) pos(3) pos(4)+0.05]); %slightly stretch y-axis
    xlabel('Distance from terminus (km)','fontsize',18); ylabel('Elevation (m)','fontsize',18); 
    subplot(subV_yr);
    set(subV_yr,'xlim',xlims,'xtick',xticks,'xticklabels',[]); 
    pos = get(subV_yr,'position'); set(subV_yr,'position',[pos(1) pos(2)-0.05 pos(3) pos(4)+0.05]); %slightly stretch y-axis
    ylabel('Speed (m/d)','fontsize',18); 
    subplot(subV_mo);
    set(subV_mo,'xlim',xlims,'xtick',xticks,'xticklabels',xticks/1000);  
    pos = get(subV_mo,'position'); set(subV_mo,'position',[pos(1) pos(2) pos(3) pos(4)+0.05]); %slightly stretch y-axis
    xlabel('Distance from terminus (km)','fontsize',18); ylabel('Speed (m/d)','fontsize',18); 
    drawnow;
    saveas(gcf,[root_dir,sitenames(j,:),'/',sitenames(j,:),'_centerline-elev-speed-terminus_subplots.png'],'png'); %save the plots
    
    %clear profile variables
    clear im im_subset LCdir zprofs zdate mean_prof seaward_idx centerline C vel_pts term_trace melmask yrs mos;
    clear DEM_num pm pt py sub* tran_* *xlims xticks *ylims yticks term_ref *_dist *_reldist *pos;
    
%     %save the structure with the centerline data
%     save([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat'],'MP','-v7.3');
%     
%     %create "climatologies" of seasonal size distributions for each site
%     D = readtable([root_dir,sitenames(j,:),'/',sitenames(j,:),'-iceberg-distribution-timeseries.csv']);
%     MP(j).D.area = D.('Area_m_2_');
%     MP(j).D.area_binwidth = D.('AreaBinwidth_m_2_');
%     %create matrices of the size distributions & their dates
%     for p = 3:size(D,2)
%         berg_dist(p-2,:) = table2array(D(:,p))';
%         berg_datestring(p-2,:) = D.Properties.VariableNames{p}(7:14);
%         berg_mo(p-2) = str2num(berg_datestring(p-2,5:6));
%     end
%     berg_dist(berg_dist==0) = NaN;
%     %create seasonal average distributions across all observations
%     %(NORMALIZED BY THE TOTAL MELANGE AREA!!)
%     dist_fig = figure; set(dist_fig,'position',[850 850 1200 600]);
%     subd1 = subplot(1,2,1); subd2 = subplot(1,2,2);
%     berg_normdist = berg_dist./sum((berg_dist.*MP(j).D.area'),2,'omitnan');
%     MP(j).D.months(1,:) = [12,1,2]; MP(j).D.months(2,:) = [3,4,5]; MP(j).D.months(3,:) = [6,7,8]; MP(j).D.months(4,:) = [9,10,11]; 
%     MP(j).D.bergs = NaN(4,size(berg_dist,2));
%     for p = 1:4
%         subplot(subd1);
%         ps(p) = loglog(NaN,NaN,'-','color',mo_cmap(p*3-2,:),'linewidth',1.5); hold on;
%         mo_ref = find(berg_mo==MP(j).D.months(p,1) | berg_mo==MP(j).D.months(p,2) | berg_mo==MP(j).D.months(p,3));
%         if ~isempty(mo_ref)
%             %plot the individual distributions colored by season
%             subplot(subd2);
%             for l = mo_ref
%                 %plot the full distribution profile for each date
%                 loglog(MP(j).D.area,berg_normdist(l,:),'-','color',mo_cmap(p*3-2,:),'linewidth',1.5); hold on;
%                 
%                 %find the first NaN in the dataset and replace all larger
%                 %size distribution counts with NaNs to avoid weird jumps in
%                 %averaged seasonal profiles due to data gaps
%                 nan_start = find(isnan(berg_normdist(l,:)) == 1,1,'first');
%                 berg_normdist(l,nan_start:end) = NaN;
%                 clear nan_start;
%             end
%             
%             %calculate the average, area-normalized distributions
%             MP(j).D.bergs(p,:) = nanmean(berg_normdist(mo_ref,:),1);
%             
%             %plot the averaged distributions
%             subplot(subd1);
%             loglog(MP(j).D.area,MP(j).D.bergs(p,:),'-','color',mo_cmap(p*3-2,:),'linewidth',2); hold on;
% 
%         end
%         clear mo_ref;
%     end
%     seas_leg = legend(ps,'DJF','MAM','JJA','SON');
%     set(subd1,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',18);
%     xlabel('Surface area (m^2)','fontsize',18); ylabel('Normalized iceberg count','fontsize',18); 
%     subplot(subd1); grid on; drawnow;
%     set(subd2,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',18);
%     xlabel('Surface area (m^2)','fontsize',18); 
%     subplot(subd2); grid on; drawnow;
%     saveas(dist_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'_seasonal-iceberg-distribution_loglog.png'],'png'); %save the plot
%     clear D berg_* ps;
%     
%     %save the structure with the seasonal distribution data
%     save([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat'],'MP','-v7.3');
%     close all;
    
end

%% create a gif that loops through and plots iceberg size distributions
%distribution for a site, then adds the fall distribution (if it exists)
%with arrows showing change between seasons for small (~100 m^2) and big
%(~10^5 m^2) icebergs

%find size indices for arrows
si = find(MP(1).D.area <= 100,1,'last');
bi = find(MP(1).D.area <= 10^5,1,'last');

%create the figure
dist_fig = figure; set(dist_fig,'position',[50 850 800 600]);
loglog(MP(1).D.area,MP(sort_ind(1)).D.bergs(2,:),'-','color',mo_cmap(2*3-2,:),'linewidth',2); hold on;
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',20); grid on; 
xlabel('Surface area (m^2)','fontsize',20); ylabel('Normalized iceberg count','fontsize',20); 
drawnow;
nimages = 1;
for j = 1:length(MP)
    if sum(~isnan(MP(sort_ind(j)).D.bergs(2,:))) > 0
        %plot spring
        loglog(MP(sort_ind(j)).D.area,MP(sort_ind(j)).D.bergs(2,:),'-','color',mo_cmap(2*3-2,:),'linewidth',2); hold on;
        title(MP(sort_ind(j)).name);
        drawnow;
        frame = getframe(dist_fig);
        gif_im{nimages} = frame2im(frame); nimages = nimages+1;
        
        for p = 3:4
            %plot summer and/or spring then fall depending on data
            if sum(~isnan(MP(sort_ind(j)).D.bergs(p,:))) > 0
                loglog(MP(sort_ind(j)).D.area,MP(sort_ind(j)).D.bergs(2,:),'-','color',mo_cmap(2*3-2,:),'linewidth',2); hold on;
                title(MP(sort_ind(j)).name);
                drawnow;
                
                %later distribution
                loglog(MP(sort_ind(j)).D.area,MP(sort_ind(j)).D.bergs(p,:),'-','color',mo_cmap(p*3-2,:),'linewidth',2); hold on;
                drawnow;
                %             frame = getframe(dist_fig);
                %             gif_im{nimages} = frame2im(frame); nimages = nimages+1;
                
                %add arrows for small iceberg abundance change
                sp1 = [MP(1).D.area(si) MP(sort_ind(j)).D.bergs(2,si)];
                sp2 = [MP(1).D.area(si) MP(sort_ind(j)).D.bergs(p,si)];
                sdp = sp2-sp1;
                quiver(sp1(1),sp1(2),sdp(1),sdp(2),0,'linewidth',2,'color','k'); hold on;
                if sdp(2) > 0
                    plot(sp2(1),sp2(2)-0.2*sp2(2),'^k','markerfacecolor','k'); hold on;
                else
                    plot(sp2(1),sp2(2)+0.2*sp2(2),'vk','markerfacecolor','k'); hold on;
                end
                clear sp1 sp2 dsp;
                %add arrows for large iceberg abundance change
                sp1 = [MP(1).D.area(bi) MP(sort_ind(j)).D.bergs(2,bi)];
                sp2 = [MP(1).D.area(bi) MP(sort_ind(j)).D.bergs(p,bi)];
                sdp = sp2-sp1;
                quiver(sp1(1),sp1(2),sdp(1),sdp(2),0,'linewidth',2,'color','k'); hold on;
                if sdp(2) > 0
                    plot(sp2(1),sp2(2)-0.2*sp2(2),'^k','markerfacecolor','k'); hold on;
                else
                    plot(sp2(1),sp2(2)+0.2*sp2(2),'vk','markerfacecolor','k'); hold on;
                end
                clear sp1 sp2 dsp;
                drawnow;
                frame = getframe(dist_fig);
                gif_im{nimages} = frame2im(frame); nimages = nimages+1;
                
                %clear the plot
                cla;
            end
        end
        cla;
        title('');
        frame = getframe(dist_fig);
        gif_im{nimages} = frame2im(frame); nimages = nimages+1;
    end
end
close;

filename = [root_dir,'Greenland-seasonal-iceberg-size-distributions_timeseries.gif']; % Specify the output file name
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
close all; clear A; clear frame gif_im map;

%plot all the seasonal distributions on a single figure
dist_fig = figure; set(dist_fig,'position',[50 50 1800 500]);
sub1 = subplot(1,3,1); sub2 = subplot(1,3,2); sub3 = subplot(1,3,3);
for j = 1:length(MP)
        %plot spring
        subplot(sub1);
        loglog(MP(sort_ind(j)).D.area,MP(sort_ind(j)).D.bergs(2,:),'-','color',mo_cmap(2*3-2,:),'linewidth',2); hold on;
        drawnow;
        
        %plot summer
        subplot(sub2);
        loglog(MP(sort_ind(j)).D.area,MP(sort_ind(j)).D.bergs(3,:),'-','color',mo_cmap(3*3-2,:),'linewidth',2); hold on;
        
        %plot fall
        subplot(sub3);
        loglog(MP(sort_ind(j)).D.area,MP(sort_ind(j)).D.bergs(4,:),'-','color',mo_cmap(4*3-2,:),'linewidth',2); hold on;
        drawnow;

end
subplot(sub1);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',20); grid on; 
xlabel('Surface area (m^2)','fontsize',20); ylabel('Normalized iceberg count','fontsize',20);
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
subplot(sub2);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'yticklabel',[],'fontsize',20); grid on; 
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
subplot(sub3);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'yticklabel',[],'fontsize',20); grid on; 
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
drawnow;
saveas(dist_fig,[root_dir,'Greenland-seasonal-iceberg-distribution_loglog.png'],'png'); %save the plot


%% create seasonally-averaged thickness & speed profiles and plot with seasonal size distributions from near the terminus & near the seaward margin
close all; drawnow;

for j = 1:length(MP)
    sitenames(j,:) = MP(j).name;
    seasons = MP(j).D.months;
    
    %set-up centerline coordinate system
    C = shaperead([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'_centerline.shp']);
    centerline_dist = [0, cumsum(sqrt((C.X(2:end)-C.X(1:end-1)).^2 + (C.Y(2:end)-C.Y(1:end-1)).^2))]';
    for l = 1:length(MP(j).V.X)
        dists = sqrt((MP(j).V.X(l)-C.X).^2 + (MP(j).V.Y(l)-C.Y).^2);
        tran_dist(l) = centerline_dist(find(dists == min(dists)));
        clear dists;
    end
    %find the reference terminus position
    for p = 1:length(MP(j).Z.date)
        term_trace(p) = MP(j).T.qualflag(p);
        zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
    end
    term_ref = find(abs(zdate-2020.66) == min(abs(zdate(term_trace==1)-2020.66))); %use terminus delineation closest to Aug. 2020 as the centerline reference
    tran_reldist = MP(j).T.centerline(1,term_ref) - tran_dist;
    centerline_reldist = MP(j).T.centerline(1,term_ref) - centerline_dist;
    seaward_idx = find(~isnan(nanmean(MP(j).Z.Zavg,2))==1,1,'first');
    inland_idx = find(~isnan(nanmean(MP(j).Z.Zavg,2))==1,1,'last');
    
    %create average 
    for p = 1:size(MP(j).Z.Zavg,2)
        datest = char(MP(j).Z.date(p));
        yrs(p) = str2num(datest(1:4)); mos(p) = str2num(datest(5:6));
        clear datest;
    end
    
    %create an average velocity profile using all velocity
    %observations with a mid-date within 15 days of the DEM
    for p = 1:length(years)
        yr_idx = find(yrs == years(p));
        if ~isempty(yr_idx)
            z_yr = MP(j).Z.Zavg(:,yr_idx); mos_yr = mos(yr_idx);
            for k = 1:4
                z_seas(:,k,p) = nanmean(z_yr(:,ismember(mos_yr,seasons(k,:))==1),2);
            end
            clear z_yr mos_yr;
        end
        clear yr_idx
    end
    clear mos yrs;
    z_seas(z_seas==0) = NaN; %remove mysterious zeros (need to look into their source!)
    
    %load the velocity timeseries for the transect-centerline intersection
    %points & create seasonally-averaged climatologies
    vel_pts = dir([root_dir,sitenames(j,:),'/velocities/']);
    for i = 1:length(vel_pts)
        if contains(vel_pts(i).name,'velocity')
            pt_ref = str2num(vel_pts(i).name(end-5:end-4));
            
            %read the file
            V = readtable([root_dir,sitenames(j,:),'/velocities/',vel_pts(i).name]);
            
            %filter out all the velocities based on temporal resolution
            short_dts = find(V.days_dt<60); %get rid of all velocities with coarse temporal resolution
            vel_dates = V.mid_date(short_dts); vel_dts = V.days_dt(short_dts);
            mos = str2num(datestr(vel_dates,'mm')); yrs = str2num(datestr(vel_dates,'yyyy'));
            vels = V.velocity_m_yr_(short_dts); vels(vels == 0) = NaN;
            
            %create an average velocity profile using all velocity
            %observations with a mid-date within 15 days of the DEM
            for p = 1:length(years)
                yr_idx = find(yrs == years(p));
                if ~isempty(yr_idx)
                    vels_yr = vels(yr_idx); mos_yr = mos(yr_idx);
                    for k = 1:4
                        vels_seas(pt_ref,k,p) = nanmean(vels_yr(ismember(mos_yr,seasons(k,:))==1));
                    end
                    clear vels_yr mos_yr;
                end
                clear yr_idx
            end
            
            clear pt_ref V short_dts vel_dates vel_dts mos yrs vels;
        end
    end
    
    %loop through the subsetted size distributions and compute seasonal
    %averages from at least the bin closest to the 2020 terminus (used as
    %the centerline reference) and the most seaward bin, adding an
    %intermediate bin for glaciers with long fjords
    Dsubs = dir([root_dir,MP(j).name,'/',MP(j).name,'*-iceberg-distribution-subsets.csv']);
    %only use complete datasets for the averaging (data was not saved for
    %every bin in each DEM!)
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(j).name,'/',Dsubs(p).name]);
        bins(p) = size(D,2)-2;
    end
    bin_no = max(bins);
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(j).name,'/',Dsubs(p).name]);
%         if (size(D,2)-2) == bin_no
            berg_mo(p) = str2num(Dsubs(p).name(9:10));
            
            %first 2 columns are area & bin width, so 3+ are data from points along the centerline
            berg_dist_seaward(p,:) = table2array(D(:,3))';
            berg_dist_inland(p,:) = table2array(D(:,size(D,2)))';
%             if size(D,2) >= 5; berg_dist_mid(p,:) = table2array(D(:,round((size(D,2)-2)/2)+2))'; end
%         else
%             berg_mo(p) = NaN;
%             berg_dist_seaward(p,:) = NaN(size(table2array(D(:,1))'));
%             berg_dist_inland(p,:) = NaN(size(table2array(D(:,1))'));
%             if size(D,2) >= 5; berg_dist_mid(p,:) = NaN(size(table2array(D(:,1))')); end
%         end
        clear D;
    end
    w = who;
    berg_dist_seaward(berg_dist_seaward==0) = NaN; berg_normdist_seaward = berg_dist_seaward./sum((berg_dist_seaward.*MP(j).D.area'),2,'omitnan');
    berg_dist_inland(berg_dist_inland==0) = NaN; berg_normdist_inland = berg_dist_inland./sum((berg_dist_inland.*MP(j).D.area'),2,'omitnan');
%     if ismember('berg_dist_mid',w)
%         berg_dist_mid(berg_dist_mid==0) = NaN; 
%         berg_normdist_mid = berg_dist_mid./sum((berg_dist_mid.*MP(j).D.area'),2,'omitnan');
%     end
    
    %create seasonal average distributions for each fjord bin: 1 = seaward,
    %2 = inland, 3 = intermediate (if it exists)
    for k = 1:4
        berg_normdist_seas(k,:,1) = nanmean(berg_normdist_seaward(ismember(berg_mo,seasons(k,:))==1,:),1);
        berg_normdist_seas(k,:,2) = nanmean(berg_normdist_inland(ismember(berg_mo,seasons(k,:))==1,:),1);
%         if ismember('berg_dist_mid',w); berg_normdist_seas(k,:,3) = nanmean(berg_normdist_mid(ismember(berg_mo,seasons(k,:))==1,:),1); end
    end
    
        
    %plot a figure of seasonally-averaged elevations on top, velocities in the middle,
    %& subsetted seasonal size distributions on the bottom
    af_fig = figure; set(af_fig,'position',[50 850 1200 600]);
    if ismember('berg_dist_mid',w)
        subz = subplot(3,3,[1:3]); subv = subplot(3,3,[4:6]); 
    else
        subz = subplot(3,2,[1:2]); subv = subplot(3,2,[3:4]); 
    end
    subplot(subz);
    for k = 1:4
        zdist = tran_reldist(seaward_idx:end);
        zmean = nanmean(z_seas(seaward_idx:end,k,:),3);
        zmax = (nanmean(z_seas(seaward_idx:end,k,:),3)+std(z_seas(seaward_idx:end,k,:),0,3,'omitnan')); zmax_idx = find(~isnan(zmax)==1);
        zmin = (nanmean(z_seas(seaward_idx:end,k,:),3)-std(z_seas(seaward_idx:end,k,:),0,3,'omitnan')); zmin_idx = find(~isnan(zmin)==1);
        if sum(~isnan(zmean)) > 0
            fill([zdist(zmax_idx), fliplr(zdist(zmin_idx))]',[zmax(zmax_idx); flipud(zmin(zmin_idx))],...
                mo_cmap(k*3-2,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
            pz(k) = plot(zdist(~isnan(zmean))',zmean(~isnan(zmean)),'-','color',mo_cmap(k*3-2,:),'linewidth',3); hold on;
        else
           pz(k) = plot(NaN,NaN,'-','color',mo_cmap(k*3-2,:),'linewidth',3); hold on; 
        end
        z_ylim(k) = max(zmax);
        clear zdist zmean zmax* zmin*;
    end
    seas_leg = legend(pz,'DJF','MAM','JJA','SON');
    set(gca,'fontsize',20); grid on; drawnow;
    set(subz,'xlim',[min(MP(j).T.centerline(1,term_ref)-MP(j).T.centerline),max(ceil(tran_reldist(seaward_idx:inland_idx)/1000)*1000)],...
        'xticklabel',[],'ylim',[0 ceil(max(z_ylim))]); ylims = get(gca,'ylim'); yticks = get(gca,'ytick');
    if length(yticks) == 2; set(gca,'ytick',[0 round(max(yticks)/2) 2*round(max(yticks)/2)],'ylim',ylims); end
    clear yticks ylims;
    ylabel('Elevation (m)','fontsize',20); %xlabel('Distance from terminus (km)','fontsize',20); 
    drawnow;
    subplot(subv);
    for k = 1:4
        %             fill([tran_reldist(seaward_idx:end), fliplr(tran_reldist(seaward_idx:end))]',...
        %                 [max(vels_seas(seaward_idx:end,k,:)./365,[],3); flipud(min(vels_seas(seaward_idx:end,k,:)./365,[],3))],...
        %                 mo_cmap(k*3-2,:),'FaceAlpha',0.25); hold on;
        vdist = tran_reldist(seaward_idx:end);
        vmean = nanmean(vels_seas(seaward_idx:end,k,:),3);
        vmax = (nanmean(vels_seas(seaward_idx:end,k,:),3)+std(vels_seas(seaward_idx:end,k,:),0,3,'omitnan'))./365; vmax_idx = find(~isnan(vmax)==1);
        vmin = (nanmean(vels_seas(seaward_idx:end,k,:),3)-std(vels_seas(seaward_idx:end,k,:),0,3,'omitnan'))./365; vmin_idx = find(~isnan(vmin)==1);
        fill([vdist(vmax_idx), fliplr(vdist(vmin_idx))]',[vmax(vmax_idx); flipud(vmin(vmin_idx))],...
            mo_cmap(k*3-2,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
        pv(k) = plot(vdist(~isnan(vmean))',vmean(~isnan(vmean))./365,'-','color',mo_cmap(k*3-2,:),'linewidth',3); hold on;
        clear vdist vmean vmax* vmin*;
    end
    set(gca,'fontsize',20); grid on; drawnow;
    set(subv,'xlim',[min(MP(j).T.centerline(1,term_ref)-MP(j).T.centerline),max(ceil(tran_reldist(seaward_idx:inland_idx)/1000)*1000)]); 
    ylims = get(gca,'ylim'); set(gca,'ylim',[0 max(ylims)]); clear ylims;
    xlabel('Distance from terminus (km)','fontsize',20); ylabel('Speed (m/d)','fontsize',20);
    pos = get(subv,'position'); set(subv,'position',[pos(1) pos(2)+0.05 pos(3) pos(4)]);
    drawnow;
    
    
    %add subplots beneath the figure that generally show the location of
    %the subsetted size distributions
%     if ismember('berg_dist_mid',w)
%         %inland-most bin
%         subplot(3,3,7);
%         for k = 1:4
%             loglog(MP(j).D.area,berg_normdist_seas(k,:,2),'-','color',mo_cmap(k*3-2,:),'linewidth',2); hold on;
%         end
%         set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],...
%             'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',20); grid on;
%         xlabel('Surface area (m^2)','fontsize',20); ylabel('Normalized iceberg count','fontsize',20);
%         text(10000,10^-1,'near-terminus','fontsize',20);
%         
%         %seaward-most bin
%         subplot(3,3,9);
%         for k = 1:4
%             loglog(MP(j).D.area,berg_normdist_seas(k,:,1),'-','color',mo_cmap(k*3-2,:),'linewidth',2); hold on;
%         end
%         set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],...
%             'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',20); grid on;
%         xlabel('Surface area (m^2)','fontsize',20); 
%         text(10000,10^-1,'seaward','fontsize',20);
%         
%         %intermediate bin
%         subplot(3,3,8);
%         for k = 1:4
%             loglog(MP(j).D.area,berg_normdist_seas(k,:,3),'-','color',mo_cmap(k*3-2,:),'linewidth',2); hold on;
%         end
%         set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],...
%             'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',20); grid on;
%         xlabel('Surface area (m^2)','fontsize',20); 
%         text(10000,10^-1,'mid-melange','fontsize',20);
%         
%     else
        %inland-most bin
        subplot(3,2,5);
        for k = 1:4
            loglog(MP(j).D.area,berg_normdist_seas(k,:,2),'-','color',mo_cmap(k*3-2,:),'linewidth',2); hold on;
        end
        set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],...
            'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',20); grid on;
        xlabel('Surface area (m^2)','fontsize',20); ylabel('Iceberg count','fontsize',20);
        text(10000,10^-1,'near-terminus','fontsize',20);
        
        %seaward-most bin
        subplot(3,2,6);
        for k = 1:4
            loglog(MP(j).D.area,berg_normdist_seas(k,:,1),'-','color',mo_cmap(k*3-2,:),'linewidth',2); hold on;
        end
        set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],...
            'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',20); grid on;
        xlabel('Surface area (m^2)','fontsize',20); 
        text(10000,10^-1,'seaward','fontsize',20);
%     end
    drawnow;
    saveas(af_fig,[root_dir,MP(j).name,'/',MP(j).name,'_seasonal-speed-size_profiles.png'],'png'); %save the plots
    
    
    %refresh
    clear berg_dist* berg_normdist* C centerline* Dsubs *idx seaward_* term_* tran_* sub* vel_* vels* w zdate berg_mo bins bin_no z_* pos pz pv;
    close all; drawnow;
end


%% download S2 images with preautorift & make gifs for Alison & Zachariae
%identfy the files (run this section for each study site, specifying
%site-specific parameters below before each rerun)

%site-specific info
site_abbrev = 'ZIM'; site_name = 'Zachariae Isstrom';
ref_image = 'S2A_27XWH_20200802_3_L2A_B08_clipped.tif'; %Alison = 'S2A_21XWC_20200731_1_L2A_B08_clipped.tif'

%load the data
load([root_dir,site_abbrev,'/',site_abbrev,'-melange-masks.mat']); %load the melange mask file
ims = dir([root_dir,site_abbrev,'/images/S2/','S*B08_clipped.tif']); im_refs = []; im_dates = [];
for k = 1:length(ims)
    if contains(ims(k).name,'_2019') || contains(ims(k).name,'_2020')
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



