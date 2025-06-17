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
years = 2011:1:2023;

%identify the site folders
cd(root_dir);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.') && length(sites(i).name) == 3
        sitenames = [sitenames; sites(i).name];
    end
end

%loop through the folders & extract info
disp('Width-averaged elevation of each transect & corresponding centerline velocity');
MP = struct;
for j = 1:length(sitenames)
    disp(sitenames(j,:));
    cd([root_dir,sitenames(j,:),'/']);
    MP(j).name = sitenames(j,:);

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
        if length(xis)>1
            MP(j).T.X(1,p) = xis(end); MP(j).T.Y(1,p) = yis(end);
            MP(j).T.centerline(1,p) = tran_dist(iis(2,end))+sqrt((MP(j).T.X(1,p)-MP(j).V.X(iis(2,end))).^2 + (MP(j).T.Y(1,p)-MP(j).V.Y(iis(2,end))).^2);
        end
        clear xis yis iis;
    end
    
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
    subZ_yr = subplot(2,3,1); yr_cmap = cmocean('matter',length(years)+2); yr_cmap = yr_cmap(2:end,:);
    subZ_mo = subplot(2,3,4); mo_cmap = cmocean('phase',12);
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
    grid on; drawnow;
    subplot(subZ_mo); 
    plot(tran_dist(end_ind:end)'-tran_dist(end_ind),mean_prof(end_ind:end),'-','color','k','linewidth',3); hold on; 
    leg_mo = legend(pm,num2str([1:1:12]')); leg_mo.Location = 'northwest';
    grid on; drawnow;
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
    grid on; drawnow;
    subplot(subV_mo); 
    plot(tran_dist(end_ind:end)'-tran_dist(end_ind),nanmean(MP(j).V.V(end_ind:end,:),2),'-','color','k','linewidth',3); hold on; 
    grid on; drawnow;
    Vxlims = get(subV_mo,'xlim');
    
    %create a terminus position timeseries
%     term_fig = figure;
    subT = subplot(2,3,[3,6]);
    for p = 1:length(MP(j).Z.date)
        mo = str2num(MP(j).Z.date{p}(5:6));
        plot(MP(j).T.centerline(p)-tran_dist(end_ind),zdate(p),'x','color',mo_cmap(mo,:),'linewidth',2); hold on;
        clear mo;
    end
    set(gca,'xlim',[0 max(MP(j).T.centerline)-tran_dist(end_ind)]); xticks = get(gca,'xtick');
    set(gca,'xtick',xticks,'xticklabels',xticks/1000);
    xlabel('Centerline distance (km)');
    %UPDATE ALL THE PLOTS SO THAT THE CENTERLINE DISTANCE IS MOVING AWAY
    %FROM THE MOST-RETREATED TERMINUS POSITION OR THE AVERAGE POSITION,
    %INSTEAD OF INLAND FROM THE OCEAN (BECAUSE THE START IS REALLY
    %ARBITRARY)
    grid on; drawnow;
%     Txlims = get(gca,'xlim');
    
    %standardize the x-limits on the plots
%     set(subZ_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subZ_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
%     set(subV_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subV_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
    set(subZ_yr,'xlim',[0 max(MP(j).T.centerline)],'xtick',xticks,'xticklabels',xticks/1000); 
    set(subZ_mo,'xlim',[0 max(MP(j).T.centerline)],'xtick',xticks,'xticklabels',xticks/1000);
    set(subV_yr,'xlim',[0 max(MP(j).T.centerline)],'xtick',xticks,'xticklabels',xticks/1000); 
    set(subV_mo,'xlim',[0 max(MP(j).T.centerline)],'xtick',xticks,'xticklabels',xticks/1000);    
    drawnow;
    
    
    %clear profile variables
    clear zprofs zdate *_cmap mean_prof end_ind C vel_pts;
    
    %compile size distributions
    
    
end




