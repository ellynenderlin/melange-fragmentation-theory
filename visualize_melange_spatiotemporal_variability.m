%%% Create plots & maps of GrIS melange data

clearvars; close all; warning off;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/inpoly2/');
addpath('/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-fragmentation-code/');

%specify paths
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-melt/';
years = 2011:1:2023;

%identify the site folders
cd(root_dir);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.')
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
    figure; set(gcf,'position',[50 50 1000 1000]);
    subZ_yr = subplot(2,2,1); yr_cmap = cmocean('matter',length(years)+2); yr_cmap = yr_cmap(2:end,:);
    subZ_mo = subplot(2,2,3); mo_cmap = cmocean('phase',12);
    mean_prof = nanmean(zprofs,2); end_ind = find(~isnan(mean_prof)==1,1,'first');
    for p = 1:size(zprofs,2)
        yr = str2num(MP(j).Z.date{p}(1:4)); mo = str2num(MP(j).Z.date{p}(5:6));
        subplot(subZ_yr);
        plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',zprofs(end_ind:end,p),'-','color',yr_cmap(yr-min(years)+1,:),'linewidth',2); hold on; 
        subplot(subZ_mo);
        plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',zprofs(end_ind:end,p),'-','color',mo_cmap(mo,:),'linewidth',2); hold on;
        clear yr mo;
    end
    subplot(subZ_yr);
    title([sitenames(j,:),' elevation profiles'])
    plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',mean_prof(end_ind:end),'-','color','k','linewidth',1.5); hold on; 
    grid on; drawnow;
    subplot(subZ_mo); 
    plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',mean_prof(end_ind:end),'-','color','k','linewidth',1.5); hold on; 
    grid on; drawnow;

    %load the shapefile of transect-centerline intersections used to
    %extract the velocity timeseries
    C = readtable([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'_centerline_2000m-interval.csv']);
    MP(j).V.X = C.Easting_m_; MP(j).V.Y = C.Northing_m_; clear C;

    %load the velocity timeseries for the transect-centerline intersection
    %points and plot a velocity profile with the closest mid-date to each
    %elevation profile (if after 2013, when Landsat 8 was launched)
    vel_pts = dir([root_dir,sitenames(j,:),'/velocities/']);
    for i = 1:length(vel_pts)
        if contains(vel_pts(i).name,'velocity')
            pt_ref = str2num(vel_pts(i).name(end-5:end-4));

            %read the file
            V = readtable([root_dir,sitenames(j,:),'/velocities/',vel_pts(i).name]);

            %filter out all the velocities with really low temporal resolution (>90 days)
            short_dts = find(V.days_dt<90);
            vel_dates = V.mid_date(short_dts); vel_dts = V.days_dt(short_dts);
            vels = V.velocity_m_yr_(short_dts);

            %convert datetime to a decimal date
            for k = 1:length(vel_dates)
                decidate(k) = convert_to_decimaldate(vel_dates(k),'yyyy-MM-dd HH:mm:ss.SSS');
            end

            %find the closest velocity mid-date to each elevation profile
            for p = 1:length(MP(j).Z.date)
                zdate = convert_to_decimaldate(char(MP(j).Z.date(p)));
                if zdate > 2013
                    datediff = abs(zdate - decidate);
                    MP(j).V.date(pt_ref,p) = string(datetime(vel_dates(find(abs(datediff) == min(datediff),1,'first')),'Format','yyyyMMdd'));
                    MP(j).V.dt(pt_ref,p) = vel_dts(find(abs(datediff) == min(datediff),1,'first'));
                    MP(j).V.V(pt_ref,p) = vels(find(abs(datediff) == min(datediff),1,'first'));
                    clear datediff;
                else
                    MP(j).V.date(pt_ref,p) = NaN; MP(j).V.dt(pt_ref,p) = NaN;
                    MP(j).V.V(pt_ref,p) = NaN;
                end
                clear zdate;
            end

            clear pt_ref V short_dts vel_dates vel_dts vels decidate;
        end
    end

    %create velocity profile plots
    for p = 1:length(MP(j).Z.date)
        subV_yr = subplot(2,2,2); subV_mo = subplot(2,2,4);
        yr = str2num(MP(j).Z.date{p}(1:4)); mo = str2num(MP(j).Z.date{p}(5:6));
        subplot(subV_yr);
        plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',MP(j).V.V(end_ind:end,p),'--','color',yr_cmap(yr-min(years)+1,:),'linewidth',2); hold on;
        subplot(subV_mo);
        plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',MP(j).V.V(end_ind:end,p),'--','color',mo_cmap(mo,:),'linewidth',2); hold on;
        clear yr mo;
    end
    subplot(subV_yr); 
    title([sitenames(j,:),' velocity profiles'])
    grid on; drawnow;
    subplot(subV_mo); grid on; drawnow;

    clear zprofs *_cmap mean_prof end_ind C vel_pts;
end




