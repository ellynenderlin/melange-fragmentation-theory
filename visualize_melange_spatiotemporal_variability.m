%%% Create plots & maps of GrIS melange data

clearvars; close all; warning off;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/inpoly2/');
% addpath('/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-fragmentation-code/');
addpath('/Users/ellynenderlin/Research/miscellaneous/melange-fragmentation-code/');

%specify root path
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';
% root_dir = '/Volumes/Jokulhaup_5T/Greenland-melange/';

%define custom parameters for size distributions
transect_inc = 1000; %distance between transects along the centerline (meters)
ARcomp.best.autoALL = 2; % iceberg aspect ratio 

%customize visualization
years = 2011:1:2023; yr_cmap = cmocean('matter',length(years)+1); yr_cmap = yr_cmap(2:end,:);
seasons = [12,1,2;3,4,5;6,7,8;9,10,11]; season_names = {'DJF','MMA','JJA','SON'};
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

%give a warning that the velocity data need to be extracted using a Google
%Colab notebook and then downloaded locally
%(Import_itslive_point_timeseries.ipynb)
disp('Extract ITS_LIVE velocities via Import_itslive_point_timeseries.ipynb & save locally before continuing!')

%% option to reload the existing data if any data have been saved
if exist([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']) == 2
    reload = questdlg('Reload the existing centerline data & start from there?',...
        'data reload','1) Yes: reload','2) No: start fresh','1) Yes: reload');
    switch reload
        case '1) Yes: reload'
            load([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']);
            if length(MP) < size(sitenames,1)
                site_start = length(MP)+1;
                SKIP = 0;
                disp(['...restarting on site #',num2str(site_start),' (',sitenames(site_start,:),')']);
            else
                disp('data reloaded and dataset is fully processed');
                replot = questdlg('Are you only (re-)plotting main figures?',...
                    'replot','1) Yes: plot only','2) No','1) Yes: plot only');
                switch replot
                    case '1) Yes: plot only'
                        % site_start = 1;
                        SKIP = 1;
                    case '2) No'
                        SKIP = 0;
                end
            end
        case '2) No: start fresh'
            site_start = 1;
            SKIP = 0;
            MP = struct;
    end
end
disp('Move on to the next subsection in order to compile data & generate figures.')

%% loop through the folders & extract info
if SKIP == 0

disp('Creating plots of elevation, velocity, and terminus position...');
for j = site_start:length(sitenames) %default: site_start:length(sitenames)
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
    
    
    %extract the transect_inc from the files (only works if only one
    %version is saved for each site but it's needed to account for adaptive
    %transect spacing depending on site size)
    shp_files = dir([root_dir,site_abbrev,'/shapefiles/',site_abbrev,'*.shp']);
    for j = 1:length(shp_files)
        if contains(shp_files(j).name,['transects_'])
            transect_inc = str2num(shp_files(j).name(end-8:end-5));
        end
    end
    %load the shapefile of transect-centerline intersections 
    C = readtable([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'-centerline_',num2str(transect_inc),'m-interval.csv'],"VariableNamingRule","preserve");
    MP(j).V.X = C.Easting_m_; MP(j).V.Y = C.Northing_m_; clear C;
    clear C;
    C = shaperead([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'-centerline.shp']);
    centerline_dist = [0, cumsum(sqrt((C.X(2:end)-C.X(1:end-1)).^2 + (C.Y(2:end)-C.Y(1:end-1)).^2))]';
    for l = 1:length(MP(j).V.X)
        dists = sqrt((MP(j).V.X(l)-C.X).^2 + (MP(j).V.Y(l)-C.Y).^2);
        tran_dist(l) = centerline_dist(find(dists == min(dists)));
        clear dists;
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

    %open the CSV of time-stamped elevation transects
    T = readtable([sitenames(j,:),'-transects_elevations.csv'],"VariableNamingRule","preserve");
    T_headers = T.Properties.VariableNames;
    datestart = strfind(string(T_headers(3)),'20');
    for p = 3:size(T,2)
        transect_date(1,p-2) = string(T_headers{p}(datestart:datestart+7));
    end

    %make sure there isn't a bad date in any of the files
    good_days = zeros(1,size(transect_date,2));
    for p = 1:size(melmask.dated,2)
        melmask_dates(p,:) = melmask.dated(p).datestring;
    end
    %first check if the melange mask file has data for dates that were not
    %used to create transects
    remove_ind = 1;
    for p = 1:size(melmask.dated,2)
        if sum(strcmp(transect_date,string(melmask_dates(p,:)))) == 0
            disp(['Removing ',melmask_dates(p,:),' from the melange mask matfile'])
            mel.uncropped = melmask.uncropped;
            mel.dated = melmask.dated(1:remove_ind-1);
            mel.dated(remove_ind:length(melmask.dated)-1) = melmask.dated(remove_ind+1:length(melmask.dated));
            clear melmask; melmask = mel; clear mel;
            save([root_dir,sitenames(j,:),'/',sitenames(j,:),'-melange-masks.mat'],'melmask','-v7.3');
        else
            remove_ind = remove_ind + 1;
        end
    end
    clear melmask_dates remove_ind;
    %then check that the transects do not have data for more dates than the
    %melange mask (which may happen if the transects are re-extracted & not all bad DEMs were deleted)
    for p = 1:size(melmask.dated,2); melmask_dates(p,:) = melmask.dated(p).datestring; end
    for p = 1:size(transect_date,2)
        good_days(1,find(strcmp(melmask_dates,transect_date(p))==1)) = 1;
    end
    remove_ind = find(good_days==0)+2; %coordinates are the first two columns so shift the bad date index by 2
    if ~isempty(remove_ind)
        T_edit = removevars(T, remove_ind);
        writetable(T_edit,[sitenames(j,:),'-transects_elevations.csv']);
        clear T; T = T_edit; clear T_edit T_headers datestart transect_date good_days remove_ind;
        %recompile good dates
        T_headers = T.Properties.VariableNames;
        datestart = strfind(string(T_headers(3)),'20');
    end
    clear transect_date remove_ind good_days;
    for p = 3:size(T,2)
        MP(j).Z.date(1,p-2) = string(T_headers{p}(datestart:datestart+7));
        datest = datetime(MP(j).Z.date{p-2},'InputFormat','yyyyMMdd');
        yrs(p-2) = year(datest); mos(p-2) = month(datest);
        clear datest;
    end
    
    %plot the image
    temp_fig = figure; set(temp_fig,'position',[50 50 1200 1200]);
    for p = 1:size(melmask.dated,2)
        imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; drawnow; hold on;
        set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)],'fontsize',20);
        plot(C.X,C.Y,'.k'); hold on;
        plot(melmask.dated(p).x,melmask.dated(p).y,'-m','linewidth',2); drawnow;
        answer = questdlg('Is the terminus delineation good?',...
            'terminus delineation','1) Yes: wiggly & good','2) No: DEM edge (straight)','3) No: wonky polygon','1) Yes: wiggly & good');
        switch answer
            case '1) Yes: wiggly & good'
                term_trace = [term_trace; 1]; melmask.dated(length(term_trace)).terminus = 1;
                MP(j).Z.termflag(p) = 1;
            case '2) No: DEM edge (straight)'
                term_trace = [term_trace; 0]; melmask.dated(length(term_trace)).terminus = 1;
                MP(j).Z.termflag(p) = 0;
            case '3) No: wonky polygon'
                removed_flag = fix_individual_melange_masks(root_dir,sitenames(j,:),melmask,melmask.dated(p).datestring);
                if strmatch(removed_flag,'removed')
                    %removed the DEM from melmask so reload it
                    load([MP(j).name,'-melange-masks.mat']);
                    DEM_num = size(melmask.dated,2); p = p-1; %reset counters to account for removed data
                else
                    term_trace = [term_trace; 1]; melmask.dated(length(term_trace)).terminus = 1;
                    MP(j).Z.termflag(p) = 1;
                end
        end
        clear answer; cla;
    end
    close(temp_fig);
    save([root_dir,sitenames(j,:),'/',sitenames(j,:),'-melange-masks.mat'],'melmask','-v7.3');

    %add terminus data from TermPicks
    % load the shapefile
    cd([root_dir,site_abbrev,'/termini/']);
    term = shaperead([site_abbrev,'_termpicks.shp']);
    for k = 1:length(term)
        term_date(k) = datenum(term(k).Date,'yyyy-mm-dd');
    end
    [~,idx] = sort(term_date);
    for k = 1:length(term)
        sorted_term(k) = term(idx(k));
    end
    clear term; term = sorted_term; clear sorted_term term_date;

    %convert format of dates
    bad_ind = [];
    for k = 1:length(term)
        %convert date to same format as used for elevation data
        if length(term(k).Month) == 1; term(k).Month = ['0',term(k).Month]; end
        if length(term(k).Day) == 1; term(k).Day = ['0',term(k).Day]; end
        YYYYMMDD(k) = string([term(k).Year,term(k).Month,term(k).Day]);
        %convert to datetimes
        try
            datestr(k) = datetime(term(k).Date,'InputFormat','yyyy-mm-dd');
        catch
            bad_ind = [bad_ind,k];
        end
    end
    %remove data with erroneous dates
    term(bad_ind) = [];
    YYYYMMDD(bad_ind) = [];

    %find the intersection of each terminus trace with the centerline and save to the structure
    for k = 1:length(term)
        [xis,yis,iis] = polyxpoly(term(k).X,term(k).Y,MP(j).V.X,MP(j).V.Y);
        MP(j).T.date(k) = YYYYMMDD(k);
        if ~isempty(xis)
            MP(j).T.termX(1,k) = xis(end); MP(j).T.termY(1,k) = yis(end);
        else
            MP(j).T.termX(1,k) = NaN; MP(j).T.termY(1,k) = NaN;
        end
        clear xis yis iis;
    end
    MP(j).T.date(isnan(MP(j).T.termX)) = []; 
    MP(j).T.termY(isnan(MP(j).T.termX)) = []; MP(j).T.termX(isnan(MP(j).T.termX)) = []; 
    %clear variables
    clear term YYYYMMDD datestr idx site_abbrev;

    %intersect the centerline with each time-stamped melange outline
    for p = 1:size(melmask.dated,2)
        zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
        [xis,yis,iis] = polyxpoly(melmask.dated(p).x,melmask.dated(p).y,MP(j).V.X,MP(j).V.Y);
        MP(j).Z.termX(1,p) = xis(end); MP(j).Z.termY(1,p) = yis(end);
        if term_trace(p) == 1 %terminus was mapped from the DEM
            if length(xis)>1
                MP(j).Z.termdist(1,p) = tran_dist(iis(end,2))+sqrt((MP(j).Z.termX(1,p)-MP(j).V.X(iis(end,2))).^2 + (MP(j).Z.termY(1,p)-MP(j).V.Y(iis(end,2))).^2);
            else
                MP(j).Z.termdist(1,p) = tran_dist(iis(1,2))+sqrt((MP(j).Z.termX(1,p)-MP(j).V.X(iis(1,2))).^2 + (MP(j).Z.termY(1,p)-MP(j).V.Y(iis(1,2))).^2);
            end
        else %terminus was cut-off in the DEM so don't record the centerline intersection
            MP(j).Z.termdist(1,p) = NaN;
        end
        clear xis yis iis;
    end
    term_ref = find(abs(zdate-2020.66) == min(abs(zdate(term_trace==1)-2020.66))); %use terminus delineation closest to Aug. 2020 as the centerline reference
    tran_reldist = MP(j).Z.termdist(1,term_ref) - tran_dist;
    centerline_reldist = MP(j).Z.termdist(1,term_ref) - centerline_dist;

    %Extract width-averaged elevations from the transects
    %find the NaNs in the coordinate pairs to identify each transect
    coords = table2array(T(:,1:2));
    nan_inds = find(isnan(coords(:,1))==1);
    
    %extract the data
    tran_ind = 1; k = 1; zprofs = [];
    while k < size(T,1)
        if k == 1
            %extract full coordinates for transect
            MP(j).Z.transectX(tran_ind,1:2) = coords([1,nan_inds(k)-1],2);
            MP(j).Z.transectY(tran_ind,1:2) = coords([1,nan_inds(k)-1],1);

            %extract the centroid coordinates
            MP(j).Z.centerX(tran_ind,1) = nanmean(table2array(T(1:nan_inds(k)-1,2)));
            MP(j).Z.centerY(tran_ind,1) = nanmean(table2array(T(1:nan_inds(k)-1,1)));

            %calculate the mean elevation for each date
            ztemp = table2array(T(1:nan_inds(k)-1,3:end)); ztemp(ztemp<3) = NaN;
            MP(j).Z.transectZavg(tran_ind,:) = nanmean(ztemp); clear ztemp;

            %add to a temp matrix for plotting
            zprofs(tran_ind,:) = MP(j).Z.transectZavg(tran_ind,:);

            tran_ind = tran_ind+1;
        elseif ismember(k,nan_inds)
            if ~ismember(k+1,nan_inds) %check for back-to-back NaNs
                %extract full coordinates for transect
                MP(j).Z.transectX(tran_ind,1:2) = coords([k+1,nan_inds(find(nan_inds==k)+1)-1],2);
                MP(j).Z.transectY(tran_ind,1:2) = coords([k+1,nan_inds(find(nan_inds==k)+1)-1],1);

                %extract the centroid coordinates
                MP(j).Z.centerX(tran_ind,1) = nanmean(table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,2)));
                MP(j).Z.centerY(tran_ind,1) = nanmean(table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,1)));

                %calculate the median elevation for each date
                ztemp = table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,3:end)); ztemp(ztemp<3) = NaN;
                MP(j).Z.transectZavg(tran_ind,:) = nanmean(ztemp); clear ztemp;

                %add to a temp matrix for plotting
                zprofs(tran_ind,:) = MP(j).Z.transectZavg(tran_ind,:);

                tran_ind = tran_ind+1;
            end
        end
        k = k+1;
    end
    clear T T_headers datestart coords nan_inds;


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
        MP(j).Z.transectZavg(MP(j).Z.transectZavg==0) = NaN; %remove mysterious zeros
        for k = 1:length(years)
            yr_idx = find(yrs == years(k));
            plot(tran_reldist(seaward_idx:end)',nanmean(MP(j).Z.transectZavg(seaward_idx:end,yr_idx),2),'-','color',yr_cmap(years(k)-min(years)+1,:),'linewidth',2); hold on; 
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
    set(gca,'fontsize',20); grid on; drawnow;
    subplot(subZ_mo); 
    plot(tran_reldist(seaward_idx:end)',mean_prof(seaward_idx:end),'-','color','k','linewidth',3); hold on; 
    pos = get(subZ_mo,'position'); 
    leg_mo = legend(pm,num2str([1:1:12]')); leg_mo.Location = 'eastoutside';
    set(subZ_mo,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]); %shift plot back to pre-legend location
    set(gca,'fontsize',20); grid on; drawnow;
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
        yr = year(datetime(MP(j).Z.date{p},'InputFormat','yyyyMMdd')); mo = month(datetime(MP(j).Z.date{p},'InputFormat','yyyyMMdd'));
%         subplot(subV_yr);
%         plot(tran_reldist(seaward_idx:end)',MP(j).V.V(seaward_idx:end,p)./365,'-','color',yr_cmap(yr-min(years)+1,:),'linewidth',2); hold on;
        subplot(subV_mo);
        plot(tran_reldist(seaward_idx:end)',MP(j).V.V(seaward_idx:end,p)./365,'-','color',mo_cmap(mo,:),'linewidth',2); hold on;
        clear yr mo;
    end
    subplot(subV_yr); 
    plot(tran_reldist(seaward_idx:end)',nanmean(MP(j).V.V(seaward_idx:end,:),2)./365,'-','color','k','linewidth',3); hold on; 
    title([sitenames(j,:),' velocity profiles'])
    set(gca,'fontsize',20); grid on; drawnow;
    subplot(subV_mo); 
    plot(tran_reldist(seaward_idx:end)',nanmean(MP(j).V.V(seaward_idx:end,:),2)./365,'-','color','k','linewidth',3); hold on; 
    set(gca,'fontsize',20); grid on; drawnow;
    Vxlims = get(subV_mo,'xlim');
    
    %create a terminus position timeseries
%     term_fig = figure;
    subT = subplot(2,3,[3,6]);
    for p = 1:length(MP(j).Z.date)
%         if term_trace(p) == 1 %don't plot terminus delineations that are the DEM edge, not the true terminus
            mo = str2num(MP(j).Z.date{p}(5:6));
            plot(tran_dist(seaward_idx)-MP(j).Z.termdist(p),zdate(p),'x','color',mo_cmap(mo,:),'linewidth',2); hold on;
            clear mo;
%         end
    end
    set(subT,'xlim',[min(tran_dist(seaward_idx)-MP(j).Z.termdist),max(ceil(tran_reldist(seaward_idx:inland_idx)/1000)*1000)],'ylim',[min(years) max(years)]); 
    xlims = get(subT,'xlim'); xticks = get(subT,'xtick');
    set(gca,'xtick',xticks,'xticklabels',xticks/1000,'fontsize',20);
    xlabel('Distance from terminus (km)','fontsize',20); ylabel('Year','fontsize',20); 
    grid on; drawnow;
%     Txlims = get(gca,'xlim');
    
    %standardize the x-limits on the plots
%     set(subZ_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subZ_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
%     set(subV_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subV_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
    subplot(subZ_yr);
    set(subZ_yr,'xlim',xlims,'xtick',xticks,'xticklabels',[]); 
    pos = get(subZ_yr,'position'); set(subZ_yr,'position',[pos(1) pos(2)-0.05 pos(3) pos(4)+0.05]); %slightly stretch y-axis
    ylabel('Elevation (m)','fontsize',20); 
    subplot(subZ_mo);
    set(subZ_mo,'xlim',xlims,'xtick',xticks,'xticklabels',xticks/1000);
    pos = get(subZ_mo,'position'); set(subZ_mo,'position',[pos(1) pos(2) pos(3) pos(4)+0.05]); %slightly stretch y-axis
    xlabel('Distance from terminus (km)','fontsize',20); ylabel('Elevation (m)','fontsize',20); 
    subplot(subV_yr);
    set(subV_yr,'xlim',xlims,'xtick',xticks,'xticklabels',[]); 
    pos = get(subV_yr,'position'); set(subV_yr,'position',[pos(1) pos(2)-0.05 pos(3) pos(4)+0.05]); %slightly stretch y-axis
    ylabel('Speed (m/d)','fontsize',20); 
    subplot(subV_mo);
    set(subV_mo,'xlim',xlims,'xtick',xticks,'xticklabels',xticks/1000);  
    pos = get(subV_mo,'position'); set(subV_mo,'position',[pos(1) pos(2) pos(3) pos(4)+0.05]); %slightly stretch y-axis
    xlabel('Distance from terminus (km)','fontsize',20); ylabel('Speed (m/d)','fontsize',20); 
    drawnow;
    saveas(gcf,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-centerline-elev-speed-terminus_subplots.png'],'png'); %save the plots
    
    %clear profile variables
    clear im im_subset LCdir zprofs zdate mean_prof seaward_idx centerline C vel_pts term_trace melmask yrs mos;
    clear DEM_num pm pt py sub* tran_* *xlims xticks *ylims yticks term_ref *_dist *_reldist *pos;
    
    %save the structure with the centerline data
    save([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat'],'MP','-v7.3');

    %create "climatologies" of seasonal size distributions for each site
    D = readtable([root_dir,sitenames(j,:),'/',sitenames(j,:),'-iceberg-distribution-timeseries.csv'],"VariableNamingRule","preserve");
    MP(j).D.area = D.('Area_m_2_');
    MP(j).D.area_binwidth = D.('AreaBinwidth_m_2_');
    %create matrices of the size distributions & their dates
    for p = 3:size(D,2)
        bergdist(p-2,:) = table2array(D(:,p))';
        berg_datestring(p-2,:) = D.Properties.VariableNames{p}(7:14);
        berg_mo(p-2) = str2num(berg_datestring(p-2,5:6));
    end
    bergdist(bergdist==0) = NaN;
    %create seasonal average distributions across all observations
    %(NORMALIZED BY THE TOTAL MELANGE AREA!!)
    dist_fig = figure; set(dist_fig,'position',[850 850 1200 600]);
    subd1 = subplot(1,2,1); subd2 = subplot(1,2,2);
    berg_normdist = bergdist./sum((bergdist.*MP(j).D.area'),2,'omitnan');
    MP(j).D.months(1,:) = seasons(1,:); MP(j).D.months(2,:) = seasons(2,:); MP(j).D.months(3,:) = seasons(3,:); MP(j).D.months(4,:) = seasons(4,:); 
    MP(j).D.bergs = NaN(4,size(bergdist,2));
    for p = 1:4
        subplot(subd1);
        ps(p) = loglog(NaN,NaN,'-','color',mo_cmap(p*3-2,:),'linewidth',1.5); hold on;
        mo_ref = find(berg_mo==MP(j).D.months(p,1) | berg_mo==MP(j).D.months(p,2) | berg_mo==MP(j).D.months(p,3));
        if ~isempty(mo_ref)
            %plot the individual distributions colored by season
            subplot(subd2);
            for l = mo_ref
                %plot the full distribution profile for each date
                loglog(MP(j).D.area,berg_normdist(l,:),'-','color',mo_cmap(p*3-2,:),'linewidth',1.5); hold on;

                %find the first NaN in the dataset and replace all larger
                %size distribution counts with NaNs to avoid weird jumps in
                %averaged seasonal profiles due to data gaps
                nan_start = find(isnan(berg_normdist(l,:)) == 1,1,'first');
                berg_normdist(l,nan_start:end) = NaN;
                clear nan_start;
            end

            %calculate the average, area-normalized distributions
            MP(j).D.bergs(p,:) = nanmean(berg_normdist(mo_ref,:),1);

            %plot the averaged distributions
            subplot(subd1);
            loglog(MP(j).D.area,MP(j).D.bergs(p,:),'-','color',mo_cmap(p*3-2,:),'linewidth',2); hold on;

        end
        clear mo_ref;
    end
    seas_leg = legend(ps,'DJF','MAM','JJA','SON');
    set(subd1,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',20);
    xlabel('Surface area (m^2)','fontsize',20); ylabel('Normalized iceberg count','fontsize',20); 
    subplot(subd1); grid on; drawnow;
    set(subd2,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',20);
    xlabel('Surface area (m^2)','fontsize',20); 
    subplot(subd2); grid on; drawnow;
    saveas(dist_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-seasonal-iceberg-distribution_loglog.png'],'png'); %save the plot
    clear D berg_* ps;

    %save the structure with the seasonal distribution data
    save([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat'],'MP','-v7.3');
    close all;
    
end

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
text(10^5.8,10^-1,'spring','fontsize',20)
xlabel('Surface area (m^2)','fontsize',20); ylabel('Normalized iceberg count','fontsize',20);
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
subplot(sub2);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'yticklabel',[],'fontsize',20); grid on;
text(10^5.8,10^-1,'summer','fontsize',20)
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
subplot(sub3);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'yticklabel',[],'fontsize',20); grid on; 
text(10^5.8,10^-1,'fall','fontsize',20)
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
drawnow;
saveas(dist_fig,[root_dir,'Greenland-seasonal-iceberg-distribution_loglog.png'],'png'); %save the plot

%% create a combined terminus position timeseries
for j = 1:length(MP)
    close all; drawnow;
    ts_fig = figure; set(ts_fig,'position',[50 50 600 1200]);

    %extract the transect_inc from the files (only works if only one
    %version is saved for each site but it's needed to account for adaptive
    %transect spacing depending on site size)
    shp_files = dir([root_dir,site_abbrev,'/shapefiles/',site_abbrev,'*.shp']);
    for j = 1:length(shp_files)
        if contains(shp_files(j).name,['transects_'])
            transect_inc = str2num(shp_files(j).name(end-8:end-5));
        end
    end

    %load the centerline & create a distance vector
    C = shaperead([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'-centerline.shp']);
    centerline_dist = [0, cumsum(sqrt((C.X(2:end)-C.X(1:end-1)).^2 + (C.Y(2:end)-C.Y(1:end-1)).^2))]';
    for l = 1:length(MP(j).V.X)
        dists = sqrt((MP(j).V.X(l)-C.X).^2 + (MP(j).V.Y(l)-C.Y).^2);
        tran_dist(l) = centerline_dist(find(dists == min(dists)));
        clear dists;
    end

    %DEM data
    for p = 1:length(MP(j).Z.date)
        term_trace(p) = MP(j).Z.termflag(p);
        zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
        zdatetime(p) = datetime(MP(j).Z.date{p},'InputFormat','yyyyMMdd');
        zyrs(p) = year(zdatetime(p)); zmos(p) = month(zdatetime(p));
    end
    mean_prof = nanmean(MP(j).Z.transectZavg,2);
    seaward_idx = find(~isnan(mean_prof)==1,1,'first');
    inland_idx = find(~isnan(mean_prof)==1,1,'last');

    %intersect the centerline with each time-stamped melange outline
    for p = 1:length(MP(j).T.date)
        Tdate(p) = convert_to_decimaldate(char(MP(j).T.date(p)));
        Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        Tyrs(p) = year(Tdatetime(p)); Tmos(p) = month(Tdatetime(p));
        dists = sqrt((MP(j).T.termX(1,p)-C.X).^2 + (MP(j).T.termY(1,p)-C.Y).^2);
        [~,term_ind] = min(dists);
        MP(j).T.termdist(1,p) = centerline_dist(term_ind);

        clear dists term_ind;
    end
    MP(j).T.termdist(Tyrs < min(years) | Tyrs > max(years)) = NaN;
    % term_ref = find(abs(zdate-2020.66) == min(abs(zdate(term_trace==1)-2020.66))); %use terminus delineation closest to Aug. 2020 as the centerline reference
    % tran_reldist = MP(j).Z.termdist(1,term_ref) - tran_dist;
    % centerline_reldist = MP(j).Z.termdist(1,term_ref) - centerline_dist;

    %compile the data for filling in terminus gaps
    term_decidates = [Tdate, zdate(term_trace ==1)]; term_dates = [Tdatetime, zdatetime(term_trace ==1)]; 
    term_X = [MP(j).T.termX, MP(j).Z.termX(term_trace ==1)];
    term_Y = [MP(j).T.termY, MP(j).Z.termY(term_trace ==1)];
    term_dists = [MP(j).T.termdist, MP(j).Z.termdist(term_trace ==1)];

    %plot the TermPicks data first, then DEM
    for p = 1:length(MP(j).T.date)
        plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).T.termdist(p),Tdate(p),'x',...
            'color',mo_cmap(Tmos(p),:),'linewidth',2,'markersize',10); hold on;
    end
    for p = 1:length(MP(j).Z.date)
        if term_trace(p) == 0
            % MP(j).Z.termX(p) = NaN; MP(j).Z.termY(p) = NaN; MP(j).Z.termdist(p) = NaN;

            dts = 365*(zdate(p)-term_decidates);
            [mindt,minind] = min(abs(dts));
            if mindt < 30
                MP(j).Z.termX(p) = term_X(minind); MP(j).Z.termY(p) = term_Y(minind);
                MP(j).Z.termdist(p) = term_dists(minind);
                plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).Z.termdist(p),zdate(p),'s',...
                    'color','k','linewidth',1,'markerfacecolor','none','markersize',10); hold on;
                disp(['Filled terminus for ',MP(j).name,' ',char(MP(j).Z.date(p)),' with ',char(term_dates(minind))]);
            else
                disp(['Need terminus data for ',MP(j).name,' ',char(MP(j).Z.date(p))]);
            end
            clear mindt minind;
        else
             plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).Z.termdist(p),zdate(p),'s',...
                 'color',mo_cmap(zmos(p),:),'linewidth',1,'markerfacecolor','none','markersize',10); hold on;
        end
       
    end
    set(gca,'xlim',[0,max([MP(j).T.termdist,MP(j).Z.termdist])-min([MP(j).T.termdist,MP(j).Z.termdist])],'ylim',[min(years) max(years)]);
    xlims = get(gca,'xlim'); xticks = get(gca,'xtick');
    set(gca,'xtick',xticks,'xticklabels',xticks/1000,'fontsize',20);
    xlabel('Distance from most-retreated terminus (km)','fontsize',20); ylabel('Year','fontsize',20);
    grid on; drawnow;
    
    %save the figure & updated terminus positions
    % save([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat'],'MP','-v7.3');
    saveas(ts_fig,[root_dir,MP(j).name,'/',MP(j).name,'-terminus-timeseries-plot.png'],'png'); %save the plots
    % uiwait %advance only after figure is closed

    clear C centerline_dist dts inland_idx seaward_idx mean_prof Tdate* term_* Tmos tran_dist Tyrs zdate* zmos zyrs;
end
close all;
disp('Done compiling terminus timeseries');


%% create seasonally-averaged thickness & speed profiles and plot with seasonal size distributions from near the terminus & near the seaward margin
close all; drawnow;

%Thickness parameters:
zcutoff = 3; %elevation threshold below which to ignore iceberge (m)
rho_i = 900; %ice density (kg/m^3)
rho_w = 1026; %water density (kg/m^3)
Hcutoff = round((rho_w/(rho_w-rho_i))*zcutoff); %H threshold for figure naming

%Velocity parameters: ideally only use velocities from time periods with
%DEMs to create the seasonal profiles but you can also use all velocities
%from that season (even if there is no DEM for that season-year combination)
vfilter = 'annual'; %options: 'annual' OR 'all'

%decide whether to sample data at (1) fixed locations or (2) relative
%locations with respect to the terminus positions in the DEMs
sampling = 'dated'; % options: 'dated' OR 'fixed';

disp(['Creating profile plots ignoring icebergs thinner than ',num2str(Hcutoff),'m & using ',vfilter,' speeds'])

%iterate
for j = 1:length(MP)
    sitenames(j,:) = MP(j).name; seasons = MP(j).D.months;
    disp(sitenames(j,:));
    
    %set-up centerline coordinate system
    C = shaperead([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'-centerline.shp']);
    centerline_dist = [0, cumsum(sqrt((C.X(2:end)-C.X(1:end-1)).^2 + (C.Y(2:end)-C.Y(1:end-1)).^2))]';
    for l = 1:length(MP(j).V.X)
        dists = sqrt((MP(j).V.X(l)-C.X).^2 + (MP(j).V.Y(l)-C.Y).^2);
        tran_dist(l) = centerline_dist(find(dists == min(dists)));
        clear dists;
    end
    %isolate the dates for data sorting
    for p = 1:length(MP(j).Z.date)
        term_trace(p) = MP(j).Z.termflag(p);
        zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
        datest = datetime(MP(j).Z.date{p},'InputFormat','yyyyMMdd');
        zyrs(p) = year(datest); zmos(p) = month(datest);
        clear datest;
    end
    % term_ref = find(abs(zdate-2020.66) == min(abs(zdate(term_trace==1)-2020.66))); %use terminus delineation closest to Aug. 2020 as the centerline reference

    %loop through subset size distributions & find the first and last
    %non-NaN columns to identify melange extent
    Dsubs = dir([root_dir,MP(j).name,'/',MP(j).name,'*-iceberg-distribution-subsets.csv']);
    size_classes = [];
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(j).name,'/',Dsubs(p).name],"VariableNamingRule","preserve");
        name_split = split(Dsubs(p).name,'-',2);
        berg_dates(p) = datetime(name_split{2},'Inputformat','yyyyMMdd'); clear name_split;

        %identify observational limits along the centerline
        berg_nos = table2array(D(:,3:end)); berg_nos(berg_nos==0) = NaN;
        size_classes(p,:) = sum(~isnan(berg_nos),1);
        seaward_ext(p) = find(size_classes(p,:)>0,1,'first');
        inland_ext(p) = find(size_classes(p,:)>0,1,'last')+1;

        %create average thickness profile
        Havg(p,:) = sum((2/ARcomp.best.autoALL)*sqrt(table2array(D(zcutoff+1:end,1))./pi()).*table2array(D(zcutoff+1:end,3:end)),1)./sum(table2array(D(zcutoff+1:end,3:end)),1);
        packing(p,:) = sum(table2array(D(zcutoff+1:end,1)).*table2array(D(zcutoff+1:end,3:end)),1)./sum(table2array(D(:,1)).*table2array(D(:,3:end)),1);
    end

    %assign seaward & inland sampling limits based on the selected sampling
    %strategy (sampling = ['fixed','dated']
    inland_meanidx = round(mean(inland_ext(term_trace==1)));
    seaward_meanidx = round(mean(seaward_ext));
    if contains(sampling,'fix')
        inland_idx = inland_meanidx.*ones(size(seaward_ext));
        seaward_idx = seaward_meanidx.*ones(size(seaward_ext)); 
    elseif contains(sampling,'date')
        %if the terminus was visible in the DEM, use it as your reference
        inland_idx = inland_ext; 
        
        %if the terminus was NOT visible in the DEM, use a combination of
        %the DEM terminus timeseries and TermPicks timeseries to
        %approximate the terminus location (gaps were filled for termini
        %with observations within 60 days in the previous section)
        for p = 1:length(term_trace)
            if term_trace(p) == 0
                if ~isnan(MP(j).Z.termdist(p))
                    dists = sqrt((MP(j).Z.termX(p)-MP(j).V.X).^2 + (MP(j).Z.termY(p)-MP(j).V.Y).^2);
                    [~,ia] = sort(dists);
                    inland_idx(p) = max(ia(1:2)); term_trace(p) = 1; %use terminus traces from imagery
                    clear dists ia;
                else
                    inland_idx(p) = NaN;
                    disp(['Need terminus data for ',MP(j).name,' ',char(MP(j).Z.date(p))])
                end
            end
        end
        seaward_idx = inland_idx-(inland_meanidx-seaward_meanidx);
    else
        error('Profile sampling strategy is not properly defined: sampling but be ''fixed'' or ''dated''')
    end
    tran_reldist = (tran_dist(inland_meanidx) - tran_dist)+transect_inc/2;
    centerline_reldist = (tran_dist(inland_meanidx) - centerline_dist)+transect_inc/2;

    %create an annual average seasonal thickness profile
    Zfilt = MP(j).Z.transectZavg; Zfilt(MP(j).Z.transectZavg==0) = NaN; Zfilt(:,term_trace==0) = NaN;
    Havg(term_trace==0,:) = NaN;
    for p = 1:length(years)
        yr_idx = find(zyrs == years(p));
        if ~isempty(yr_idx)
            mos_yr = zmos(yr_idx);

            %create a matrix of elevations relative to the terminus position
            z_profiles = NaN(max(inland_idx),length(yr_idx));
            H_profiles = NaN(length(yr_idx),max(inland_idx)-1);
            pack_profiles = NaN(length(yr_idx),max(inland_idx)-1);
            for k = 1:length(yr_idx)
                if ~isnan(inland_idx(yr_idx(k)))
                    z_temp = flipud(Zfilt(1:inland_idx(yr_idx(k)),yr_idx(k)));
                    z_profiles(1:length(z_temp),k) = z_temp; 
                    H_temp = fliplr(Havg(yr_idx(k),1:inland_idx(yr_idx(k))-1));
                    H_profiles(k,1:length(H_temp)) = H_temp; 
                    pack_temp = fliplr(packing(yr_idx(k),1:inland_idx(yr_idx(k))-1));
                    pack_profiles(k,1:length(pack_temp)) = pack_temp; 
                    clear z_temp H_temp pack_temp;
                end
            end

            %calculate seasonal averages
            for k = 1:4
                z_seas(:,k,p) = nanmean(z_profiles(:,ismember(mos_yr,seasons(k,:))==1),2);
                % H_seas(:,k,p) = (rho_w./(rho_w-rho_i))*nanmean(z_profiles(:,ismember(mos_yr,seasons(k,:))==1),2);
                H_seas(:,k,p) = nanmean(H_profiles(ismember(mos_yr,seasons(k,:))==1,:),1)';
                if ~isempty(inland_idx(ismember(mos_yr,seasons(k,:))==1))
                    MP(j).B.ref(1,k,p) = mean(inland_idx(ismember(mos_yr,seasons(k,:))==1));
                    MP(j).B.Ho(zcutoff+1,k,p) = nanmean(H_profiles(ismember(mos_yr,seasons(k,:))==1,1),1)';
                    MP(j).B.packing(zcutoff+1,k,p) = nanmean(pack_profiles(ismember(mos_yr,seasons(k,:))==1,1),1)';
                else
                    MP(j).B.ref(1,k,p) = NaN;
                    MP(j).B.Ho(zcutoff+1,k,p) = NaN;
                    MP(j).B.packing(zcutoff+1,k,p) = NaN;
                end
            end

            clear z_profiles H_profiles pack_profiles mos_yr;
        else
            %no DEM for that year so fill in all thickness-based datasets
            %used to calculate buttressing with NaNs
            MP(j).B.ref(1,1:4,p) = NaN;
            MP(j).B.Ho(zcutoff+1,1:4,p) = NaN;
            MP(j).B.packing(zcutoff+1,1:4,p) = NaN;
        end
        clear yr_idx
    end
    z_seas(z_seas==0) = NaN; H_seas(H_seas==0) = NaN; 
    zdist = 0:transect_inc:(max(inland_idx)-1)*transect_inc; 
    Hdist = zdist+transect_inc/2;

    %load the velocity timeseries for the transect-centerline intersection
    %points & create seasonally-averaged climatologies for each year
    vel_pts = dir([root_dir,sitenames(j,:),'/velocities/']);
    vel_seas = NaN(max(inland_idx),4,length(years));
    for i = 1:length(vel_pts)
        if contains(vel_pts(i).name,'velocity')
            pt_ref = str2num(vel_pts(i).name(end-5:end-4));

            if pt_ref <= max(inland_idx)

                %read the file
                V = readtable([root_dir,sitenames(j,:),'/velocities/',vel_pts(i).name]);

                %filter out all the velocities based on temporal resolution
                short_dts = find(V.days_dt<60); %get rid of all velocities with coarse temporal resolution
                vel_dates = V.mid_date(short_dts); vel_dts = V.days_dt(short_dts);
                vmos = month(vel_dates); vyrs = year(vel_dates);
                vels = V.velocity_m_yr_(short_dts); vels(vels == 0) = NaN;

                %calculate seasonal average speeds for each year at the point
                if contains(vfilter,'all')
                    %grab velocities for all years
                    for p = 1:length(years)
                        yr_idx = find(vyrs == years(p));
                        if ~isempty(yr_idx)

                            %isolate the velocities for that year
                            for k = 1:length(yr_idx)
                                v_temp(k,1) = vels(yr_idx(k));
                            end

                            %calculate seasonal average
                            mos_yr = vmos(yr_idx);
                            for k = 1:4
                                if ~isempty(zmos(zyrs==years(p)))
                                    %use the DEMs from that year to come up with the position wrt the terminus
                                    rel_ref = round(median(inland_idx(zyrs==years(p))))-pt_ref+1; %edit if monthly terminus positions are available
                                    if rel_ref >=1
                                        vel_seas(rel_ref,k,p) = nanmean(v_temp(ismember(mos_yr,seasons(k,:))==1));
                                    end
                                end
                            end

                            clear vels_yr mos_yr v_temp;
                        end
                        clear yr_idx
                    end

                elseif contains(vfilter,'annual')
                    %only use velocities from years with DEMs
                    for p = 1:length(years)
                        if sum(ismember(zyrs,years(p))) > 0 %only extract velocities from years with DEMs
                            yr_idx = find(vyrs == years(p));
                            if ~isempty(yr_idx)

                                %isolate the velocities for that year
                                for k = 1:length(yr_idx)
                                    v_temp(k,1) = vels(yr_idx(k));
                                end

                                %calculate seasonal average if there is a DEM
                                %for that year and within that season
                                mos_yr = vmos(yr_idx);
                                for k = 1:4
                                    if ~isempty(zmos(zyrs==years(p) & ismember(zmos,seasons(k,:))))
                                        %use the DEMs from that season to come
                                        %up with the position wrt the terminus
                                        rel_ref = round(median(inland_idx(zyrs==years(p) & ismember(zmos,seasons(k,:)))))-pt_ref+1; %edit if monthly terminus positions are available
                                        if rel_ref >=1
                                            vel_seas(rel_ref,k,p) = nanmean(v_temp(ismember(mos_yr,seasons(k,:))==1));
                                        end
                                    end
                                end

                                clear vels_yr mos_yr v_temp;
                            end
                            clear yr_idx
                        end
                    end
                else
                    error('Unkown velocity sampling strategy! vfilter must be ''annual'' or ''all''')
                end

                clear V short_dts vel_dates vel_dts vmos vyrs vels rel_ref;
            end
            clear pt_ref;
        end
    end
    vdist = 0:transect_inc:(max(inland_idx)-1)*transect_inc; 

    %loop through the subsetted size distributions and compute seasonal
    %averages nearest the terminus & a fixed relative distance down-fjord
    Dsubs = dir([root_dir,MP(j).name,'/',MP(j).name,'*-iceberg-distribution-subsets.csv']);
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(j).name,'/',Dsubs(p).name],"VariableNamingRule","preserve");
        yrs(p) = year(berg_dates(p)); mos(p) = month(berg_dates(p));

        %first 2 columns are area & bin width, so 3+ are data from points along the centerline
        if ~isnan(inland_idx(p))
            %full size distributions at two points
            bergdist_inland(p,:) = table2array(D(:,(inland_idx(p)-1)+2))'; %near-terminus
            bergdist_seaward_setdx(p,:) = table2array(D(:,seaward_ext(p)+2))'; %standard distance from near-terminus
        else
            bergdist_inland(p,:) = NaN(size(table2array(D(:,1))))'; %near-terminus
            bergdist_seaward_setdx(p,:) = NaN(size(table2array(D(:,1))))'; %standard distance from near-terminus
        end
        bergdist_seaward_end(p,:) = table2array(D(:,seaward_ext(p)+2))'; %seaward melange limit in DEM
        clear D;
    end
    bergdist_inland(bergdist_inland==0) = NaN; bergdist_inland_norm = bergdist_inland./sum((bergdist_inland.*MP(j).D.area'),2,'omitnan');
    bergdist_seaward_setdx(bergdist_seaward_setdx==0) = NaN; bergdist_seaward_setdx_norm = bergdist_seaward_setdx./sum((bergdist_seaward_setdx.*MP(j).D.area'),2,'omitnan');
    bergdist_seaward_end(bergdist_seaward_end==0) = NaN; bergdist_seaward_end_norm = bergdist_seaward_end./sum((bergdist_seaward_end.*MP(j).D.area'),2,'omitnan');

    %create seasonal average distributions for each fjord bin
    for k = 1:4
        bergdist_seas(k,:,1) = nanmean(bergdist_seaward_setdx_norm(ismember(mos,seasons(k,:))==1,:),1); %could change to bergdist_seaward_end_norm for moving end
        bergdist_seas(k,:,2) = nanmean(bergdist_inland_norm(ismember(mos,seasons(k,:))==1,:),1);
    end
        
    %plot a figure of seasonally-averaged elevations on top, velocities in the middle,
    %& subsetted seasonal size distributions on the bottom
    af_fig = figure; set(af_fig,'position',[50 850 1200 600]);
    subz = subplot(3,2,[1:2]); subv = subplot(3,2,[3:4]); 
    subplot(subz);
    for k = 1:4
        %relative to moving terminus
        Hmean = nanmean(H_seas(:,k,:),3);
        Hmax = (nanmean(H_seas(:,k,:),3)+std(H_seas(:,k,:),0,3,'omitnan')); 
        Hmin = (nanmean(H_seas(:,k,:),3)-std(H_seas(:,k,:),0,3,'omitnan')); 

        %plot
        Hmax_idx = find(~isnan(Hmax)==1); Hmin_idx = find(~isnan(Hmin)==1);
        if sum(~isnan(Hmean)) > 0
            % fill([zdist(Hmax_idx), fliplr(zdist(Hmin_idx))]',[Hmax(Hmax_idx); flipud(Hmin(Hmin_idx))],...
            %     mo_cmap(k*3-2,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
            % pz(k) = plot(zdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',mo_cmap(k*3-2,:),'linewidth',3); hold on;
            fill([Hdist(Hmax_idx), fliplr(Hdist(Hmin_idx))]',[Hmax(Hmax_idx); flipud(Hmin(Hmin_idx))],...
                mo_cmap(k*3-2,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
            pz(k) = plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',mo_cmap(k*3-2,:),'linewidth',3); hold on;
        else
           pz(k) = plot(NaN,NaN,'-','color',mo_cmap(k*3-2,:),'linewidth',3); hold on; 
        end
        % H_ylim(k,:) = [min(Hmin), max(Hmax)];
        H_ylim(k,:) = [min(Hmean), max(Hmean)];
        clear Hmean Hmax* Hmin*;
        % clear zdist;
    end
    seas_leg = legend(pz,season_names);
    set(gca,'fontsize',20); grid on; drawnow;
    set(subz,'xlim',[0,vdist(max([find(sum(sum(~isnan(H_seas),2),3)>0,1,'last'),find(sum(sum(~isnan(vel_seas),2),3)>0,1,'last')]))],...
        'xticklabel',[],'ylim',[floor(min(H_ylim(:,1))/10)*10 ceil(max(H_ylim(:,2))/10)*10]); 
    ylims = get(gca,'ylim'); 
    if range(ylims) <= 50
        set(gca,'ytick',[floor(min(H_ylim(:,1))/10)*10:10:ceil(max(H_ylim(:,2))/10)*10]);
    else
        set(gca,'ytick',[floor(min(H_ylim(:,1))/10)*10:20:ceil(max(H_ylim(:,2))/10)*10]);
    end
    % set(subz,'xlim',[0,floor(tran_reldist(seaward_idx)/1000)*1000+500],'xticklabel',[],'ylim',[0 ceil(max(H_ylim)/50)*50]); 
    clear ylims;
    ylabel('Thickness (m)','fontsize',20); %xlabel('Distance from terminus (km)','fontsize',20); 
    subz_pos = get(subz,'position'); 
    seas_leg.Location = 'northoutside'; seas_leg.Orientation = 'horizontal'; %move the legend
    set(subz,'position',subz_pos); %resize the subplot
    drawnow;
    subplot(subv);
    for k = 1:4
        %relative to moving terminus
        vmean = nanmean(vel_seas(:,k,:),3)./365;
        vmax = (nanmean(vel_seas(:,k,:),3)+std(vel_seas(:,k,:),0,3,'omitnan'))./365; 
        vmin = (nanmean(vel_seas(:,k,:),3)-std(vel_seas(:,k,:),0,3,'omitnan'))./365; 

        %plot
        vmax_idx = find(~isnan(vmax)==1); vmin_idx = find(~isnan(vmin)==1);
        if sum(~isnan(vmean)) ~= 0
            fill([vdist(vmax_idx), fliplr(vdist(vmin_idx))]',[vmax(vmax_idx); flipud(vmin(vmin_idx))],...
                mo_cmap(k*3-2,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
            pv(k) = plot(vdist(~isnan(vmean))',vmean(~isnan(vmean)),'-','color',mo_cmap(k*3-2,:),'linewidth',3); hold on;
        end
        clear vmean vmax* vmin*;
        % clear vdist;
    end
    set(gca,'fontsize',20); grid on; drawnow;
    set(subv,'xlim',[0,vdist(max([find(sum(sum(~isnan(H_seas),2),3)>0,1,'last'),find(sum(sum(~isnan(vel_seas),2),3)>0,1,'last')]))]);
    xticks = get(subv,'xtick'); set(subv,'xticklabel',xticks/1000); clear xticks;
    % set(subv,'xlim',[0,floor(tran_reldist(seaward_idx)/1000)*1000+500]); 
    ylims = get(subv,'ylim'); set(subv,'ylim',[0 max(ylims)]); 
    yticks = get(gca,'ytick'); 
    if length(yticks) < 3
        if range(ylims) <= 20
            set(gca,'ytick',[0:5:max(ylims)]);
        else
            set(gca,'ytick',[0:10:max(ylims)]);
        end
    end
    clear ylims yticks;
    xlabel('Distance from terminus (km)','fontsize',20); ylabel('Speed (m/d)','fontsize',20);
    pos = get(subv,'position'); set(subv,'position',[pos(1) pos(2)+0.05 pos(3) pos(4)]);
    drawnow;
    
    %add subplots beneath the figure that generally show the location of the size distributions
    %inland bin
    subplot(3,2,5);
    for k = 1:4
        loglog(MP(j).D.area,bergdist_seas(k,:,2),'-','color',mo_cmap(k*3-2,:),'linewidth',2); hold on;
    end
    set(gca,'ylim',[10^-10 10^-1],'xlim',[10^1 10^7],...
        'ytick',[10^-9,10^-5,10^-1],...
        'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',20); grid on;
    xlabel('Surface area (m^2)','fontsize',20); ylabel('Iceberg count','fontsize',20);
    text(10000,10^-2,'near-terminus','fontsize',20);

    %seaward bin
    subplot(3,2,6);
    for k = 1:4
        loglog(MP(j).D.area,bergdist_seas(k,:,1),'-','color',mo_cmap(k*3-2,:),'linewidth',2); hold on;
    end
    set(gca,'ylim',[10^-10 10^-1],'xlim',[10^1 10^7],...
        'ytick',[10^-9,10^-5,10^-1],...
        'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',20); grid on;
    xlabel('Surface area (m^2)','fontsize',20);
    text(10000,10^-2,'seaward','fontsize',20);
    drawnow;

    
    %estimate buttressing
    if contains(vfilter,'annual')
        %calculate buttressing
        for p = 1:length(years)
            for k = 1:4
                % add seasonal buttressing info to the structure
                MP(j).B.dVdx(1,k,p) = (diff(vel_seas(1:2,k,p))./diff(vdist(1:2)))/365;
                press = 0.5*rho_i*(1-(rho_i/rho_w))*9.81*MP(j).B.Ho(zcutoff+1,k,p);
                MP(j).B.butt_Meng(zcutoff+1,k,p) = press*MP(j).B.packing(zcutoff+1,k,p)*MP(j).B.Ho(zcutoff+1,k,p);
                MP(j).B.butt_Amundson(zcutoff+1,k,p) = (-2*(MP(j).B.Ho(zcutoff+1,k,p)*press*MP(j).B.dVdx(1,k,p))/((MP(j).B.dVdx(1,k,p)/0.3)+MP(j).B.dVdx(1,k,p)))+press*MP(j).B.Ho(zcutoff+1,k,p);
                clear press;
            end
        end

        %display buttressing data
        disp(' Buttressing estimated using near-terminus observations:');
        disp('   Meng Eqn (thickness- & packing-based): x10^6 N/m');
        disp(['    spring = ',num2str(round(nanmedian(MP(j).B.butt_Meng(zcutoff+1,2,:))/10^6,1)),', summer = ',num2str(round(nanmedian(MP(j).B.butt_Meng(zcutoff+1,3,:))/10^6,1)),', fall = ',num2str(round(nanmedian(MP(j).B.butt_Meng(zcutoff+1,4,:))/10^6,1))]);
        disp('   Amundson Eqn (thickness- & strainrate-based): x10^6 N/m');
        disp(['    spring = ',num2str(round(nanmedian(MP(j).B.butt_Amundson(zcutoff+1,2,:))/10^6,1)),', summer = ',num2str(round(nanmedian(MP(j).B.butt_Amundson(zcutoff+1,3,:))/10^6,1)),', fall = ',num2str(round(nanmedian(MP(j).B.butt_Amundson(zcutoff+1,4,:))/10^6,1))]);
        disp(' Thickness: m')
        disp(['    spring = ',num2str(round(nanmedian(MP(j).B.Ho(zcutoff+1,2,:)),1)),', summer = ',num2str(round(nanmedian(MP(j).B.Ho(zcutoff+1,3,:)),1)),', fall = ',num2str(round(nanmedian(MP(j).B.Ho(zcutoff+1,4,:)),1))]);
        disp(' Packing density: fractional area')
        disp(['    spring = ',num2str(round(nanmedian(MP(j).B.packing(zcutoff+1,2,:)),2)),', summer = ',num2str(round(nanmedian(MP(j).B.packing(zcutoff+1,3,:)),2)),', fall = ',num2str(round(nanmedian(MP(j).B.packing(zcutoff+1,4,:)),2))]);
        disp(' Strain rate characteristics: (m/d)/m')
        disp(['    spring = ',num2str(round(nanmedian(MP(j).B.dVdx(1,2,:)),4)),', summer = ',num2str(round(nanmedian(MP(j).B.dVdx(1,3,:)),4)),', fall = ',num2str(round(nanmedian(MP(j).B.dVdx(1,4,:)),4))]);
        disp(' ');
    end
    

    %save the data and the figure
    save([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat'],'MP','-v7.3');
    saveas(af_fig,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-speed-size_',num2str(Hcutoff),'m-Hthreshold_',vfilter,'-speeds_',sampling,'-profiles.png'],'png'); %save the plots
    % uiwait %advance only after figure is closed
    
    %refresh
    clear berg_* bergdist* berg_normdist* C centerline* D Dsubs *idx seaward_* inland_* term_* tran_* size_classes sub* Zfilt H_* Havg vel_* vels* w zdate berg_mo bins bin_no z_* pos pz pv *dist *yrs *mos seas_leg packing;
    close all; drawnow;
end


%% create overview maps for each site
cd(root_dir);
load([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']);

for j = 1:length(MP)
    close all; drawnow;
    disp(sitenames(j,:)); output_dir = [root_dir,sitenames(j,:),'/']; site_abbrev = MP(j).name;
    
    %load the melange masks
    cd([root_dir,sitenames(j,:)]);
    load([MP(j).name,'-melange-masks.mat']); %load the melange mask file

    %load the Landsat image
    LCdir = dir([root_dir,sitenames(j,:),'/LC*']); im_dir = [LCdir(1).folder,'/',LCdir(1).name,'/']; %Landsat 8 or 9 unzipped image directory for mapping
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

    %isolate terminus positions from the dated melange masks
    for p = 1:length(melmask.dated)
        [~,on] = inpolygon(melmask.dated(p).x,melmask.dated(p).y,melmask.uncropped.x,melmask.uncropped.y);
        melmask.dated(p).x(on) = []; melmask.dated(p).y(on) = [];
    end
    %grab info about the DEMs & corresponding terminus positions
    for p = 1:length(MP(j).Z.date)
        term_trace(p) = MP(j).Z.termflag(p);
        zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
        datest = char(MP(j).Z.date(p));
        zyrs(p) = str2num(datest(1:4)); zmos(p) = str2num(datest(5:6));
        for l = 1:4
            if ismember(zmos(p),seasons(l,:))
                zseas(p) = l;
            end
        end
        clear datest;
    end

    %add image-based terminus delineations
    term = shaperead([root_dir,site_abbrev,'/termini/',site_abbrev,'_termpicks.shp']);
    for k = 1:length(term)
        term_date(k) = datenum(term(k).Date,'yyyy-mm-dd');
    end
    [~,idx] = sort(term_date);
    for k = 1:length(term)
        sorted_term(k) = term(idx(k));
    end
    clear term; term = sorted_term; clear sorted_term term_date;

    %convert format of dates
    bad_ind = [];
    for k = 1:length(term)
        %convert date to same format as used for elevation data
        for l = 1:4
            if ismember(str2num(term(k).Month),seasons(l,:))
                Tseas(k) = l;
            end
        end
        if length(term(k).Month) == 1; term(k).Month = ['0',term(k).Month]; end
        if length(term(k).Day) == 1; term(k).Day = ['0',term(k).Day]; end
        YYYYMMDD(k) = string([term(k).Year,term(k).Month,term(k).Day]); Tyrs(k) = str2num(term(k).Year);
        %convert to datetimes
        try
            datestr(k) = datetime(term(k).Date,'InputFormat','yyyy-mm-dd');
        catch
            bad_ind = [bad_ind,k];
        end
    end
    %remove data with erroneous dates
    term(bad_ind) = []; 
    YYYYMMDD(bad_ind) = []; Tseas(bad_ind) = []; Tyrs(bad_ind) = []; 
    T_inds = find(Tyrs < min(zyrs) | Tyrs > max(zyrs));
    term(T_inds) = []; YYYYMMDD(T_inds) = []; Tseas(T_inds) = []; Tyrs(T_inds) = [];

    %identify the centerline points from which data were extracted and
    %create a colormap for the bins
    Dsubs = dir([root_dir,MP(j).name,'/',MP(j).name,'*-iceberg-distribution-subsets.csv']);
    size_classes = [];
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(j).name,'/',Dsubs(p).name],"VariableNamingRule","preserve");
        name_split = split(Dsubs(p).name,'-',2);
        berg_dates(p) = datetime(name_split{2},'Inputformat','yyyyMMdd'); clear name_split;

        %identify observational limits along the centerline
        berg_nos = table2array(D(:,3:end)); berg_nos(berg_nos==0) = NaN;
        size_classes(p,:) = sum(~isnan(berg_nos),1);
        seaward_ext(p) = find(size_classes(p,:)>0,1,'first');
        inland_ext(p) = find(size_classes(p,:)>0,1,'last');
        clear D;
    end
    clear Dsubs;

    %create a colormap for illustrating melange bins
    % tran_cmap = cmocean('deep',(max(inland_ext)-min(seaward_ext)+1));
    yellow_end = [255,255,178]/255; red_end = [189,0,38]/255; colorgrad = (yellow_end-red_end)/(max(inland_ext)-min(seaward_ext)+1);
    for k = 1:(max(inland_ext)-min(seaward_ext)+1)
        tran_cmap(k,:) = yellow_end-(k-1)*colorgrad;
    end
    tran_cmap = flipud(tran_cmap);

    %create an overview map
    map_fig = figure; set(map_fig,'position',[850 50 800 600]);
    imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; drawnow; hold on;
    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]);
    %plot dummy lines to create a legend
    first_term = find(term_trace == 1,1,'first');
    % for l = 1:size(yr_cmap,1)
    %     pt(l) = plot(melmask.dated(first_term).x,melmask.dated(first_term).y,...
    %         '-','color',yr_cmap(l,:),'linewidth',1.5); hold on;
    %     leg_labels(l) = cellstr(num2str(years(l)));
    % end
    for l = 1:4
        pt(l) = plot(melmask.dated(first_term).x,melmask.dated(first_term).y,...
            '-','color',mo_cmap(3*l-2,:),'linewidth',1.5); hold on;
    end
    %plot faint colored boxes for the DEM binning
    for k = min(seaward_ext):max(inland_ext)
        fill([MP(j).Z.transectX(k,1),MP(j).Z.transectX(k,2),MP(j).Z.transectX(k+1,2),MP(j).Z.transectX(k+1,1),MP(j).Z.transectX(k,1)],...
            [MP(j).Z.transectY(k,1),MP(j).Z.transectY(k,2),MP(j).Z.transectY(k+1,2),MP(j).Z.transectY(k+1,1),MP(j).Z.transectY(k,1)],...
            tran_cmap(k-min(seaward_ext)+1,:),'FaceAlpha',0.5,'EdgeColor','none'); hold on;
    end
    %plot the terminus positions
    for p = 1:size(melmask.dated,2)
        if term_trace(p) == 1
            % plot(melmask.dated(p).x,melmask.dated(p).y,'-','color',...
            %     yr_cmap(zyrs(p)-min(years)+1,:),'linewidth',1.5); hold on;
            plot(melmask.dated(p).x,melmask.dated(p).y,'-','color',...
                mo_cmap(3*zseas(p)-2,:),'linewidth',1.5); hold on;
            drawnow;
        end
    end
    for k = 1:length(term)
        in = inpolygon(term(k).X,term(k).Y,melmask.uncropped.x,melmask.uncropped.y);
        % plot(MP(j).T.X(in),MP(j).T.Y(in),'-','color',...
        %     yr_cmap(Tyrs(k)-min(years)+1,:),'linewidth',1.5); hold on;
        plot(term(k).X(in),term(k).Y(in),'-','color',...
            mo_cmap(3*Tseas(k)-2,:),'linewidth',1.5); hold on;
        clear in; drawnow;
    end
    %add centerline & transect locations
    plot(MP(j).V.X(min(seaward_ext):max(inland_ext)+1),MP(j).V.Y(min(seaward_ext):max(inland_ext)+1),'.-k','markersize',16,'linewidth',2); hold on;
    xlims = [min(MP(j).Z.transectX(min(seaward_ext):max(inland_ext)+1,:),[],'all'), max(MP(j).Z.transectX(min(seaward_ext):max(inland_ext)+1,:),[],'all')];
    ylims = [min(MP(j).Z.transectY(min(seaward_ext):max(inland_ext)+1,:),[],'all'), max(MP(j).Z.transectY(min(seaward_ext):max(inland_ext)+1,:),[],'all')];
    set(gca,'xlim',xlims,'ylim',ylims);
    xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
    set(gca,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',20);
    xlabel('Easting (km)','fontsize',20); ylabel('Northing (km)','fontsize',20);
    pos = get(gca,'position');
    if range(xlims) > 1.05*range(ylims) %short and fat map so plot the legend below
        % map_leg = legend(pt,leg_labels,'Location','southoutside',...
        %     'Orientation','horizontal','NumColumns',7);
        map_leg = legend(pt,season_names,'Location','southoutside',...
            'Orientation','horizontal');
        leg_pos = get(map_leg,'position');
        set(map_leg,'position',[leg_pos(1)+(0.5-mean([leg_pos(1) leg_pos(1)+leg_pos(3)])) 0.075 leg_pos(3) leg_pos(4)]);
        if range(xlims)/range(ylims) < 1.5
            set(gca,'position',[pos(1) 0.25 0.9*pos(3) 0.9*pos(4)]); drawnow;
        else
            set(gca,'position',[pos(1) 0.2 pos(3) pos(4)]); drawnow;
        end
    else %tall and thin map so plot the legend on the side
        % map_leg = legend(pt,leg_labels,'Location','eastoutside',...
        %     'Orientation','vertical','NumColumns',1);
        map_leg = legend(pt,season_names,'Location','eastoutside',...
            'Orientation','vertical');
        leg_pos = get(map_leg,'position');
        if leg_pos(1)+0.05 > pos(1)-0.05+pos(3)
            set(map_leg,'position',[leg_pos(1)+0.05 leg_pos(2) leg_pos(3) leg_pos(4)]);
        else
            set(map_leg,'position',[pos(1)-0.05+pos(3) leg_pos(2) leg_pos(3) leg_pos(4)]);
        end
        set(gca,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]); drawnow;
    end
    saveas(map_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-site-map.png'],'png'); %save the image
    uiwait %advance only after figure is closed

    clear xlims ylims pt;
    clear berg_* im* inland_ext LCdir leg* map* melmask zmos on pos ref* seaward_ext size_classes sort_ind T_inds term_* Tyrs Tmos *lims zyrs zdate zseas Tseas tran_cmap colorgrad term;
end

%% load TermPicks timeseries for each glacier
% cd(root_dir);
% load([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']);
% 
% for j = 1:length(MP)
%     site_abbrev = MP(j).name;
% 
%     % load the shapefile
%     cd([root_dir,site_abbrev,'/termini/']);
%     term = shaperead([site_abbrev,'_termpicks.shp']);
%     for k = 1:length(term)
%         term_date(k) = datenum(term(k).Date,'yyyy-mm-dd');
%     end
%     [~,idx] = sort(term_date);
%     for k = 1:length(term)
%         sorted_term(k) = term(idx(k));
%     end
%     clear term; term = sorted_term; clear sorted_term term_date;
% 
%     %convert format of dates
%     bad_ind = [];
%     for k = 1:length(term)
%         %convert date to same format as used for elevation data
%         if length(term(k).Month) == 1; term(k).Month = ['0',term(k).Month]; end
%         if length(term(k).Day) == 1; term(k).Day = ['0',term(k).Day]; end
%         YYYYMMDD(k) = string([term(k).Year,term(k).Month,term(k).Day]);
%         %convert to datetimes
%         try
%             datestr(k) = datetime(term(k).Date,'InputFormat','yyyy-mm-dd');
%         catch
%             bad_ind = [bad_ind,k];
%         end
%     end
%     %remove data with erroneous dates
%     term(bad_ind) = [];
%     YYYYMMDD(bad_ind) = [];
% 
% 
%     %find the intersection of each terminus trace with the centerline and
%     %save to the structure
%     for k = 1:length(term)
%         [xis,yis,iis] = polyxpoly(term(k).X,term(k).Y,MP(j).V.X,MP(j).V.Y);
%         MP(j).T.date(k) = YYYYMMDD(k);
%         if ~isempty(xis)
%             MP(j).T.termX(1,k) = xis(end); MP(j).T.termY(1,k) = yis(end);
%         else
%             MP(j).T.termX(1,k) = NaN; MP(j).T.termY(1,k) = NaN;
%         end
%         clear xis yis iis;
%     end
%     MP(j).T.date(isnan(MP(j).T.termX)) = []; 
%     MP(j).T.termY(isnan(MP(j).T.termX)) = []; MP(j).T.termX(isnan(MP(j).T.termX)) = []; 
% 
%     %save to the structure
%     save([root_dir,'GrIS-melange_centerline-elev-speed-terminus.mat'],'MP','-v7.3');
% 
%     %clear variables
%     clear term YYYYMMDD datestr idx site_abbrev;
% end



%% fit iceberg size distribution curves to average of normalized, binned, size distributions





