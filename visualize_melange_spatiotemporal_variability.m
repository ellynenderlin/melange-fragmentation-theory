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
years = 2011:1:2023; yr_cmap = cmocean('matter',length(years)+2); yr_cmap = yr_cmap(2:end,:);
mo_cmap = cmocean('phase',12); close all;

%identify the site folders
cd(root_dir);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.') && length(sites(i).name) == 3
        sitenames = [sitenames; sites(i).name];
    end
end

%option to reload the existing data if any data have been saved
if exist([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']) == 2
    reload = questdlg('Reload the existing centerline data & start from there?',...
        'data reload','1) Yes: reload','2) No: start fresh','1) Yes: reload');
    switch reload
        case '1) Yes: reload'
            load([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']);
            site_start = length(MP)+1;
        case '2) No: start fresh'
            site_start = 1;
    end
end

%loop through the folders & extract info
disp('Creating plots of elevation, velocity, and terminus position...');
MP = struct;
for j = site_start:length(sitenames)
    disp(sitenames(j,:)); output_dir = [root_dir,sitenames(j,:),'/'];
    MP(j).name = sitenames(j,:);
    
    %navigate to the study site directory
    cd(output_dir);
    LCdir = dir([root_dir,sitenames(j,:),'/LC*']); im_dir = [LCdir(1).folder,'/',LCdir(1).name,'/']; %Landsat 8 or 9 unzipped image directory for mapping
    
    %open the CSV of time-stamped elevation transects
    T = readtable([sitenames(j,:),'_transects_elevations.csv']);
    T_headers = T.Properties.VariableNames;
    datestart = strfind(string(T_headers(3)),'20');
    for p = 3:size(T,2)
        MP(j).Z.date(1,p-2) = string(T_headers{p}(datestart:datestart+7));
    end
    
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
    tran_dist = [0; cumsum(sqrt((MP(j).V.X(2:end)-MP(j).V.X(1:end-1)).^2 + (MP(j).V.Y(2:end)-MP(j).V.Y(1:end-1)).^2))];
    clear C;
    C = shaperead([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'_centerline.shp']);
    centerline_dist = [0, cumsum(sqrt((C.X(2:end)-C.X(1:end-1)).^2 + (C.Y(2:end)-C.Y(1:end-1)).^2))]';
    
    %intersect the centerline with each time-stamped melange outline
    for p = 1:size(melmask.dated,2)
        zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
        [xis,yis,iis] = polyxpoly(melmask.dated(p).x,melmask.dated(p).y,MP(j).V.X,MP(j).V.Y);
        MP(j).T.X(1,p) = xis(end); MP(j).T.Y(1,p) = yis(end);
        if length(xis)>1
            MP(j).T.centerline(1,p) = tran_dist(iis(end,2))+sqrt((MP(j).T.X(1,p)-MP(j).V.X(iis(end,2))).^2 + (MP(j).T.Y(1,p)-MP(j).V.Y(iis(end,2))).^2);
        else
            MP(j).T.centerline(1,p) = tran_dist(iis(1,2))+sqrt((MP(j).T.X(1,p)-MP(j).V.X(iis(1,2))).^2 + (MP(j).T.Y(1,p)-MP(j).V.Y(iis(1,2))).^2);
        end
        clear xis yis iis;
    end
    
    %identify which terminus traces are actually the edge of the DEM
    %(straight lines) then make an overview map that excludes them
    ims = dir([im_dir,'L*.TIF']);
    for k = 1:length(ims)
        if contains(ims(k).name,'B8')
            ref_image = [ims(k).folder,'/',ims(k).name];
        end
    end
    clear im_dir ims;
    %load the panchromatic Landsat scene
    [I,R] = readgeoraster(ref_image);
    im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
    im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
    im.z = double(I);
    clear I R;
    %crop the image to adjust brightnesses appropriately
    xlims = [find(im.x<=min(melmask.uncropped.x),1,'last'), find(im.x>=max(melmask.uncropped.x),1,'first')];
    ylims = [find(im.y>=max(melmask.uncropped.y),1,'last'), find(im.y<=min(melmask.uncropped.y),1,'first')];
    im_subset = im.z(min(ylims):max(ylims),min(xlims):max(xlims));
    im_subset = im_subset./max(max(im_subset));
    %plot the image
    temp_fig = figure; set(temp_fig,'position',[50 50 1200 1200]);
    for p = 1:DEM_num
        imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; drawnow; hold on;
        set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)],'fontsize',14);
        plot(C.X,C.Y,'.k'); hold on;
        plot(melmask.dated(p).x,melmask.dated(p).y,'-m','linewidth',2); drawnow;
        answer = questdlg('Is the terminus delineation good?',...
            'terminus delineation','1) Yes: wiggly & good','2) No: DEM edge (straight)','3) No: wonky polygon','1) Yes: wiggly & good');
        switch answer
            case '1) Yes: wiggly & good'
                term_trace = [term_trace; 1]; melmask.dated(length(term_trace)).terminus = 1;
            case '2) No: DEM edge (straight)'
                term_trace = [term_trace; 0]; melmask.dated(length(term_trace)).terminus = 1;
            case '3) No: wonky polygon'
                removed_flag = fix_individual_melange_masks(root_dir,sitenames(j,:),melmask,p);
                if strmatch(removed_flag,'removed')
                    %removed the DEM from melmask so reload it
                    load([MP(j).name,'-melange-masks.mat']);
                    DEM_num = size(melmask.dated,2);
                else
                    term_trace = [term_trace; 1]; melmask.dated(length(term_trace)).terminus = 1;
                end
        end
        clear answer; cla;
    end
    close(temp_fig);
    save([root_dir,sitenames(j,:),'/',sitenames(j,:),'-melange-masks.mat'],'melmask','-v7.3');
    
    %create an overview map of all the dated melange masks
    map_fig = figure; set(map_fig,'position',[650 50 600 600]);
    imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; drawnow; hold on;
    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]);
    %plot dummy lines to create a legend
    first_term = find(term_trace == 1,1,'first');
    for l = 1:size(yr_cmap,1)
        pt(l) = plot(melmask.dated(first_term).x,melmask.dated(first_term).y,...
            '-','color',yr_cmap(l,:),'linewidth',1.5); hold on;
    end
    %plot the real terminus traces
    for p = 1:size(melmask.dated,2)
        if term_trace(p) == 1
            plot(melmask.dated(p).x,melmask.dated(p).y,'-','color',...
                yr_cmap(str2num(melmask.dated(p).datestring(1:4))-min(years)+1,:),'linewidth',1.5); hold on;
            drawnow;
        end
    end
    plot(C.X,C.Y,'.k'); hold on;
    plot(melmask.uncropped.x,melmask.uncropped.y,'-b','linewidth',2); hold on;
    xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
    set(gca,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
    xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
    if range(xlims) > range(ylims) %short and fat map so plot the legend below
        map_leg = legend(pt,num2str(years'),'Location','southoutside',...
            'Orientation','horizontal','NumColumns',7);
    else %tall and thin map so plot the legend on the side
        map_leg = legend(pt,num2str(years'),'Location','eastoutside',...
            'Orientation','vertical','NumColumns',1);
    end
    drawnow;
    saveas(map_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'_site-map.png'],'png'); %save the image
    
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

            %calculate the median elevation for each date
            MP(j).Z.Zavg(tran_ind,:) = nanmean(table2array(T(1:nan_inds(k)-1,3:end)));

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
                MP(j).Z.Zavg(tran_ind,:) = nanmean(table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,3:end)));

                %add to a temp matrix for plotting
                zprofs(tran_ind,:) = MP(j).Z.Zavg(tran_ind,:);

                tran_ind = tran_ind+1;
            end
        end
        k = k+1;
    end
    clear T T_headers datestart coords nan_inds;

    %plot time-series of the width-averaged elevation profiles
    figure; set(gcf,'position',[50 50 1200 600]);
    subZ_yr = subplot(2,3,1); subZ_mo = subplot(2,3,4); 
    mean_prof = nanmean(zprofs,2); end_ind = find(~isnan(mean_prof)==1,1,'first');
    for p = 1:size(zprofs,2)
        %add dummy lines for the legend
        if p == 1
            for l = 1:size(yr_cmap,1)
                subplot(subZ_yr);
                py(l) = plot(tran_dist(end_ind:end)'-tran_dist(end_ind),zprofs(end_ind:end,p),'-','color',yr_cmap(l,:),'linewidth',2); hold on;
            end
            for l = 1:size(mo_cmap,1)
                subplot(subZ_mo);
                pm(l) = plot(tran_dist(end_ind:end)'-tran_dist(end_ind),zprofs(end_ind:end,p),'-','color',mo_cmap(l,:),'linewidth',2); hold on;
            end
        end
        
        %plot the data
        yr = str2num(MP(j).Z.date{p}(1:4)); mo = str2num(MP(j).Z.date{p}(5:6));
        subplot(subZ_yr);
        plot(tran_dist(end_ind:end)'-tran_dist(end_ind),zprofs(end_ind:end,p),'-','color',yr_cmap(yr-min(years)+1,:),'linewidth',2); hold on; 
        subplot(subZ_mo);
        plot(tran_dist(end_ind:end)'-tran_dist(end_ind),zprofs(end_ind:end,p),'-','color',mo_cmap(mo,:),'linewidth',2); hold on;
        clear yr mo;
    end
    subplot(subZ_yr);
    title([sitenames(j,:),' elevation profiles'])
    plot(tran_dist(end_ind:end)'-tran_dist(end_ind),mean_prof(end_ind:end),'-','color','k','linewidth',3); hold on; 
    leg_yr = legend(py,num2str(years')); leg_yr.Location = 'northwest';
    set(gca,'fontsize',14); grid on; drawnow;
    subplot(subZ_mo); 
    plot(tran_dist(end_ind:end)'-tran_dist(end_ind),mean_prof(end_ind:end),'-','color','k','linewidth',3); hold on; 
    leg_mo = legend(pm,num2str([1:1:12]')); leg_mo.Location = 'northwest';
    set(gca,'fontsize',14); grid on; drawnow;
    Zxlims = get(subZ_mo,'xlim');

    %load the velocity timeseries for the transect-centerline intersection
    %points and plot a velocity profile with the closest mid-date to each
    %elevation profile (if after 2013, when Landsat 8 was launched)
    vel_pts = dir([root_dir,sitenames(j,:),'/velocities/']);
    for i = 1:length(vel_pts)
        if contains(vel_pts(i).name,'velocity')
            pt_ref = str2num(vel_pts(i).name(end-5:end-4));

            %read the file
            V = readtable([root_dir,sitenames(j,:),'/velocities/',vel_pts(i).name]);

            %filter out all the velocities based on temporal resolution
            short_dts = find(V.days_dt<60); %get rid of all velocities with coarse temporal resolution
            vel_dates = V.mid_date(short_dts); vel_dts = V.days_dt(short_dts);
            vels = V.velocity_m_yr_(short_dts); vels(vels == 0) = NaN;

            %convert datetime to a decimal date
            for k = 1:length(vel_dates)
                decidate(k) = convert_to_decimaldate(vel_dates(k),'yyyy-MM-dd HH:mm:ss.SSS');
            end

            %create an average velocity profile using all velocity
            %observations with a mid-date within 15 days of the DEM
            for p = 1:length(MP(j).Z.date)
                if zdate(p) > 2013
                    datediff = abs(zdate(p) - decidate);
                    MP(j).V.date(pt_ref,p) = string(nanmean(datetime(vel_dates(find(abs(datediff) <= (15/365))),'Format','yyyyMMdd')));
                    MP(j).V.dt(pt_ref,p) = nanmean(vel_dts(find(abs(datediff) <= (15/365))));
                    MP(j).V.V(pt_ref,p) = nanmean(vels(find(abs(datediff) <= (15/365))));
                    clear datediff;
                else
                    MP(j).V.date(pt_ref,p) = NaN; MP(j).V.dt(pt_ref,p) = NaN;
                    MP(j).V.V(pt_ref,p) = NaN;
                end
            end

            clear pt_ref V short_dts vel_dates vel_dts vels decidate;
        end
    end

    %create velocity profile plots
    subV_yr = subplot(2,3,2); subV_mo = subplot(2,3,5);
    for p = 1:length(MP(j).Z.date)
        yr = str2num(MP(j).Z.date{p}(1:4)); mo = str2num(MP(j).Z.date{p}(5:6));
        subplot(subV_yr);
        plot(tran_dist(end_ind:end)'-tran_dist(end_ind),MP(j).V.V(end_ind:end,p),'-','color',yr_cmap(yr-min(years)+1,:),'linewidth',2); hold on;
        subplot(subV_mo);
        plot(tran_dist(end_ind:end)'-tran_dist(end_ind),MP(j).V.V(end_ind:end,p),'-','color',mo_cmap(mo,:),'linewidth',2); hold on;
        clear yr mo;
    end
    subplot(subV_yr); 
    plot(tran_dist(end_ind:end)'-tran_dist(end_ind),nanmean(MP(j).V.V(end_ind:end,:),2),'-','color','k','linewidth',3); hold on; 
    title([sitenames(j,:),' velocity profiles'])
    set(gca,'fontsize',14); grid on; drawnow;
    subplot(subV_mo); 
    plot(tran_dist(end_ind:end)'-tran_dist(end_ind),nanmean(MP(j).V.V(end_ind:end,:),2),'-','color','k','linewidth',3); hold on; 
    set(gca,'fontsize',14); grid on; drawnow;
    Vxlims = get(subV_mo,'xlim');
    
    %create a terminus position timeseries
%     term_fig = figure;
    subT = subplot(2,3,[3,6]);
    for p = 1:length(MP(j).Z.date)
        if term_trace(p) == 1 %don't plot terminus delineations that are the DEM edge, not the true terminus
            mo = str2num(MP(j).Z.date{p}(5:6));
            plot(MP(j).T.centerline(p)-tran_dist(end_ind),zdate(p),'x','color',mo_cmap(mo,:),'linewidth',2); hold on;
            clear mo;
        end
    end
    set(gca,'xlim',[0 max(MP(j).T.centerline)-tran_dist(end_ind)],'ylim',[min(years) max(years)]); 
    xticks = get(gca,'xtick');
    set(gca,'xtick',xticks,'xticklabels',xticks/1000,'fontsize',14);
    xlabel('Centerline distance (km)','fontsize',14); ylabel('Year','fontsize',14); 
    %UPDATE ALL THE PLOTS SO THAT THE CENTERLINE DISTANCE IS MOVING AWAY
    %FROM THE MOST-RETREATED TERMINUS POSITION OR THE AVERAGE POSITION,
    %INSTEAD OF INLAND FROM THE OCEAN (BECAUSE THE START IS REALLY
    %ARBITRARY)
    grid on; drawnow;
%     Txlims = get(gca,'xlim');
    
    %standardize the x-limits on the plots
%     set(subZ_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subZ_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
%     set(subV_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subV_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
    subplot(subZ_yr);
    set(subZ_yr,'xlim',[0 max(MP(j).T.centerline)],'xtick',xticks,'xticklabels',xticks/1000); 
    ylabel('Elevation (m)','fontsize',14); 
    subplot(subZ_mo);
    set(subZ_mo,'xlim',[0 max(MP(j).T.centerline)],'xtick',xticks,'xticklabels',xticks/1000);
    xlabel('Centerline distance (km)','fontsize',14); ylabel('Elevation (m)','fontsize',14); 
    subplot(subV_yr);
    set(subV_yr,'xlim',[0 max(MP(j).T.centerline)],'xtick',xticks,'xticklabels',xticks/1000); 
    ylabel('Speed (m/yr)','fontsize',14); 
    subplot(subV_mo);
    set(subV_mo,'xlim',[0 max(MP(j).T.centerline)],'xtick',xticks,'xticklabels',xticks/1000);  
    xlabel('Centerline distance (km)','fontsize',14); ylabel('Speed (m/yr)','fontsize',14); 
    drawnow;
    saveas(gcf,[root_dir,sitenames(j,:),'/',sitenames(j,:),'_centerline-elev-speed-terminus_subplots.png'],'png'); %save the image
    
    %clear profile variables
    clear im im_subset LCdir zprofs zdate mean_prof end_ind centerline C vel_pts term_trace melmask;
    clear DEM_num pm pt py sub* tran_* *xlims xticks *ylims yticks;
    
    %save the structure with the centerline data
    save([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat'],'MP','-v7.3');
    
    %compile size distributions
    close all;
    
end




