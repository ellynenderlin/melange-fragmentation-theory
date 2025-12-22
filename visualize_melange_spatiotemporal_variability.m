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
zthresh = 3; %cutoff elevation (m)

%customize visualization
years = 2011:1:2023; yr_cmap = cmocean('matter',length(years)+1); yr_cmap = yr_cmap(2:end,:);
seasons = [12,1,2;3,4,5;6,7,8;9,10,11]; season_names = {'DJF','MAM','JJA','SON'};
mo_cmap = cmocean('phase',12); 
seas_cmap = [5,113,176; 202,0,32; 244,165,130; 146,197,222]/255; %blue/orange/red
% seas_cmap = [194,165,207; 123,50,148; 0,136,55; 166,219,160]/255; %purple/green
close all;

%identify the site folders
cd(root_dir);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.') && length(sites(i).name) == 3
        sitenames = [sitenames; sites(i).name];
    end
end
%if you have a preferred order for the sites in plots, specify it here
geo_order = [{'ULS'},{'KOG'},{'ASG'},{'ILS'},{'UNS'},{'SAS'},{'UMS'},...
    {'KGS'},{'SEK'},{'HLG'},{'MGG'},{'MDG'},{'DJG'},{'ZIM'}];
for j = 1:length(geo_order)
    geo_ind(j) = find(contains(string(sitenames),geo_order(j)));
end
big3 = [{'SEK'},{'HLG'},{'ZIM'}];

%give a warning that the velocity data need to be extracted using a Google
%Colab notebook and then downloaded locally
%(Import_itslive_point_timeseries.ipynb)
disp('Extract ITS_LIVE velocities via Import_itslive_point_timeseries.ipynb & save locally before continuing!')

%% option to reload the existing data if any data have been saved
if exist([root_dir,'GrIS-melange-characteristics.mat']) == 2
    reload = questdlg('Reload the existing centerline data & start from there?',...
        'data reload','1) Yes: reload','2) No: start fresh','1) Yes: reload');
    switch reload
        case '1) Yes: reload'
            load([root_dir,'GrIS-melange-characteristics.mat']);
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
                        SKIP = 0; site_start = 1;
                end
            end
        case '2) No: start fresh'
            site_start = 1;
            SKIP = 0;
            MP = struct;
    end
else
    site_start = 1;
    SKIP = 0;
    MP = struct;
end
disp('Move on to the next subsection in order to compile data & generate figures.')

%% loop through the folders & extract info
if SKIP == 0
    % %redoing everything other than flagging terminus positions... load old data
    % cd(root_dir);
    % load GrIS-melange_centerline-elev-speed-terminus.mat;
    % for j = 1:length(MP)
    %     t(j).flag = MP(j).Z.termflag;
    % end
    % MP = rmfield(MP,{'Z','V'}); 
    % for j = 1:length(MP)
    %     MP(j).Z.termflag = t(j).flag;
    % end
    % clear t;

    disp('Compiling elevation, velocity, and terminus position...');
    for j = site_start:length(sitenames) %default: site_start:length(sitenames)
        disp(sitenames(j,:)); output_dir = [root_dir,sitenames(j,:),'/'];
        site_abbrev = sitenames(j,:);
        MP(j).name = sitenames(j,:); 
        % term_trace = MP(j).Z.termflag;
        
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
        for k = 1:length(shp_files)
            if contains(shp_files(k).name,['transects_'])
                transect_inc = str2num(shp_files(k).name(end-8:end-5));
            end
        end
        %load the shapefile of transect-centerline intersections
        C = readtable([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'-centerline_',num2str(transect_inc),'m-interval.csv']);
        MP(j).V.X = C.Easting_m_; MP(j).V.Y = C.Northing_m_; clear C;
        clear C;
        C = shaperead([root_dir,sitenames(j,:),'/shapefiles/',sitenames(j,:),'-centerline.shp']);
        centerline_dist = [0, cumsum(sqrt((C.X(2:end)-C.X(1:end-1)).^2 + (C.Y(2:end)-C.Y(1:end-1)).^2))]';
        for l = 1:length(MP(j).V.X)
            dists = sqrt((MP(j).V.X(l)-C.X).^2 + (MP(j).V.Y(l)-C.Y).^2);
            tran_dist(l) = centerline_dist(find(dists == min(dists)));
            clear dists;
        end

        %load the reference satellite image
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
        tebergAfig = figure; set(tebergAfig,'position',[50 50 1200 1200]);
        imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; drawnow; hold on;
        set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)],'fontsize',16);
        for p = 1:size(melmask.dated,2)
            zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
            zdatetime(p) = datetime(MP(j).Z.date{p},'InputFormat','yyyyMMdd');
            zyrs(p) = year(zdatetime(p)); zmos(p) = month(zdatetime(p));

            %plot the melange outline and flag the inland boundary as the
            %terminus if not already completed or dataset is updated
            if isfield(MP(j).Z,'termflag')
                if length(MP(j).Z.termflag) ~= length(melmask.dated)

                    %plot the melange mask for the given date
                    imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; drawnow; hold on;
                    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)],'fontsize',16);
                    plot(C.X,C.Y,'.k'); hold on;
                    plot(melmask.dated(p).x,melmask.dated(p).y,'-m','linewidth',2); drawnow;

                    %check if DEM-based delineations have already been labeled as the
                    %real terminus or the DEM edge
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
                else
                    term_trace(p) = MP(j).Z.termflag(p);
                    %plot the melange mask for the given date
                    plot(melmask.dated(p).x,melmask.dated(p).y,'-m','linewidth',2); drawnow;
                end
            end
        end
        close(tebergAfig);
        % save([root_dir,sitenames(j,:),'/',sitenames(j,:),'-melange-masks.mat'],'melmask','-v7.3');

        %add manual delineations from images
        cd([root_dir,site_abbrev,'/termini/']);
        % termfiles = dir('*term*.shp');
        termfile_names = {[site_abbrev,'_termpicks.shp'],[site_abbrev,'-termini-EPSG3413_2019-2023.shp']}; %specify consistent data sources
        for l = 1:length(termfile_names)
            % load the shapefile
            % if ismember(termfiles(l).name,termfile_names)
            term = shaperead(char(termfile_names(l)));
            for k = 1:length(term)
                term_date(k) = datenum(term(k).Date,'yyyy-mm-dd');
            end
            [~,idx] = sort(term_date);
            for k = 1:length(term)
                sorted_term(k) = term(idx(k));
            end
            clear term; term = sorted_term; clear sorted_term term_date;
            %convert format of TermPicks dates
            bad_ind = [];
            for k = 1:length(term)
                %isolate the year, moth, and day from the date
                if ~isfield(term(k),'Year')
                    term(k).Year = num2str(year(term(k).Date));
                    term(k).Month = num2str(month(term(k).Date));
                    term(k).Day = num2str(day(term(k).Date));
                else
                    if k == 1 %wipe out existing data in case format is inconsistent with what I want
                        term = rmfield(term,{'Year','Month','Day'});
                    end
                    term(k).Year = num2str(year(term(k).Date));
                    term(k).Month = num2str(month(term(k).Date));
                    term(k).Day = num2str(day(term(k).Date));
                end

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

            %find the intersection of each terminus trace with the
            %centerline and save to the structure
            if l == 1
                start_ref = 0;
                MP(j).T = rmfield(MP(j).T,{'date','termX','termY'});
                % MP(j).T.date = []; MP(j).T.termX = []; MP(j).T.termY = []; MP(j).T.termdist = [];
            else
                start_ref = length(MP(j).T.date);
            end
            for k = 1:length(term)
                [xis,yis,iis] = polyxpoly(term(k).X,term(k).Y,MP(j).V.X,MP(j).V.Y);
                MP(j).T.date(start_ref+k) = {char(YYYYMMDD(k))};
                if ~isempty(xis)
                    MP(j).T.termX(1,start_ref+k) = xis(end); MP(j).T.termY(1,start_ref+k) = yis(end);
                else
                    MP(j).T.termX(1,start_ref+k) = NaN; MP(j).T.termY(1,start_ref+k) = NaN;
                end
                clear xis yis iis;
            end
            MP(j).T.date(isnan(MP(j).T.termX)) = [];
            MP(j).T.termY(isnan(MP(j).T.termX)) = []; MP(j).T.termX(isnan(MP(j).T.termX)) = [];
            %clear variables
            clear term YYYYMMDD datestr idx;
            % end
        end

        %intersect the centerline with each time-stamped terminus position
        for p = 1:length(MP(j).T.date)
            try
                Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
            catch
                MP(j).T.date{p} = string(str2num(MP(j).T.date{p})-1);
                Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
            end
            Tdate(p) = convert_to_decimaldate(char(MP(j).T.date{p}));
            Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
            Tyrs(p) = year(Tdatetime(p)); Tmos(p) = month(Tdatetime(p));
            dists = sqrt((MP(j).T.termX(1,p)-C.X).^2 + (MP(j).T.termY(1,p)-C.Y).^2);
            [~,term_ind] = min(dists);
            MP(j).T.termdist(1,p) = centerline_dist(term_ind);

            clear dists term_ind;
        end

        %intersect the centerline with each time-stamped melange outline
        for p = 1:size(melmask.dated,2)
            [xis,yis,iis] = polyxpoly(melmask.dated(p).x,melmask.dated(p).y,MP(j).V.X,MP(j).V.Y);
            MP(j).Z.termX(1,p) = xis(end); MP(j).Z.termY(1,p) = yis(end);
            if term_trace(p) == 1 %terminus was mapped from the DEM
                if length(xis)>1
                    MP(j).Z.termdist(1,p) = tran_dist(iis(end,2))+sqrt((MP(j).Z.termX(1,p)-MP(j).V.X(iis(end,2))).^2 + (MP(j).Z.termY(1,p)-MP(j).V.Y(iis(end,2))).^2);
                else
                    MP(j).Z.termdist(1,p) = tran_dist(iis(1,2))+sqrt((MP(j).Z.termX(1,p)-MP(j).V.X(iis(1,2))).^2 + (MP(j).Z.termY(1,p)-MP(j).V.Y(iis(1,2))).^2);
                end
            else %terminus was cut-off in the DEM so didn't record the centerline intersection
                MP(j).Z.termdist(1,p) = NaN;
            end
            clear xis yis iis;
        end
        %compile the data for filling in terminus gaps
        term_decidates = [Tdate, zdate(term_trace ==1)]; term_dates = [Tdatetime, zdatetime(term_trace ==1)];
        term_X = [MP(j).T.termX, MP(j).Z.termX(term_trace == 1)];
        term_Y = [MP(j).T.termY, MP(j).Z.termY(term_trace == 1)];
        term_dists = [MP(j).T.termdist, MP(j).Z.termdist(term_trace ==1)];

        %create relative distance vectors
        for p = 1:length(term_trace)
            if term_trace(p) == 0
                %estimate the terminus position from image-based delineations
                dts = 365*(zdate(p)-term_decidates);
                [mindt,minind] = min(abs(dts));
                if mindt < 30
                    MP(j).Z.termX(p) = term_X(minind); MP(j).Z.termY(p) = term_Y(minind);
                    MP(j).Z.termdist(p) = term_dists(minind);
                    disp(['Filled terminus for ',MP(j).name,' ',char(MP(j).Z.date(p)),' with ',char(term_dates(minind))]);
                else
                    disp(['Need terminus data for ',MP(j).name,' ',char(MP(j).Z.date(p))]);
                    MP(j).Z.termdist(1,p) = NaN;
                end
                clear mindt minind;
            end
        end
        % term_ref = find(abs(zdate-2020.66) == min(abs(zdate(term_trace==1)-2020.66))); %use terminus delineation closest to Aug. 2020 as the centerline reference


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
                ztemp = table2array(T(1:nan_inds(k)-1,3:end)); ztemp(ztemp<zthresh) = NaN;
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
                    ztemp = table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,3:end)); ztemp(ztemp<zthresh) = NaN;
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

        %establish a coordinate system relative to the inland-most melange
        %elevation observations
        mean_prof = nanmean(zprofs,2);
        seaward_idx = find(~isnan(mean_prof)==1,1,'first');
        inland_idx = find(~isnan(mean_prof)==1,1,'last');
        tran_reldist = tran_dist(inland_idx) - tran_dist;
        centerline_reldist = tran_dist(inland_idx) - centerline_dist;

        %plot time-series of the width-averaged elevation profiles
        figure; set(gcf,'position',[50 50 1200 700]);
        subZ_yr = subplot(2,3,1); subZ_mo = subplot(2,3,4);
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
        set(gca,'fontsize',16); grid on; drawnow;
        subplot(subZ_mo);
        plot(tran_reldist(seaward_idx:end)',mean_prof(seaward_idx:end),'-','color','k','linewidth',3); hold on;
        pos = get(subZ_mo,'position');
        leg_mo = legend(pm,num2str([1:1:12]')); leg_mo.Location = 'eastoutside';
        set(subZ_mo,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]); %shift plot back to pre-legend location
        set(gca,'fontsize',16); grid on; drawnow;
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
        set(gca,'fontsize',16); grid on; drawnow;
        subplot(subV_mo);
        plot(tran_reldist(seaward_idx:end)',nanmean(MP(j).V.V(seaward_idx:end,:),2)./365,'-','color','k','linewidth',3); hold on;
        set(gca,'fontsize',16); grid on; drawnow;
        Vxlims = get(subV_mo,'xlim');

        %create a terminus position timeseries
        %     term_fig = figure;
        subT = subplot(2,3,[3,6]);
        for p = 1:length(MP(j).Z.date)
            %         if term_trace(p) == 1 %don't plot terminus delineations that are the DEM edge, not the true terminus
            mo = str2num(MP(j).Z.date{p}(5:6));
            plot(tran_dist(inland_idx) - MP(j).Z.termdist(p),zdate(p),'x','color',mo_cmap(mo,:),'linewidth',2); hold on;
            clear mo;
            %         end
        end
        set(subT,'xlim',[min(tran_dist(inland_idx) - MP(j).Z.termdist),max(ceil(tran_reldist(seaward_idx:inland_idx)/1000)*1000)],'ylim',[min(years) max(years)]);
        xlims = get(subT,'xlim'); xticks = get(subT,'xtick');
        set(gca,'xtick',xticks,'xticklabels',xticks/1000,'fontsize',16);
        xlabel('Distance from terminus (km)','fontsize',16); ylabel('Year','fontsize',16);
        grid on; drawnow;
        %     Txlims = get(gca,'xlim');

        %standardize the x-limits on the plots
        %     set(subZ_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subZ_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
        %     set(subV_yr,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]); set(subV_mo,'xlim',[0 min([Zxlims(2),Vxlims(2),Txlims(2)])]);
        subplot(subZ_yr);
        set(subZ_yr,'xlim',xlims,'xtick',xticks,'xticklabels',[]);
        pos = get(subZ_yr,'position'); set(subZ_yr,'position',[pos(1) pos(2)-0.05 pos(3) pos(4)+0.05]); %slightly stretch y-axis
        ylabel('Elevation (m)','fontsize',16);
        subplot(subZ_mo);
        set(subZ_mo,'xlim',xlims,'xtick',xticks,'xticklabels',xticks/1000);
        pos = get(subZ_mo,'position'); set(subZ_mo,'position',[pos(1) pos(2) pos(3) pos(4)+0.05]); %slightly stretch y-axis
        xlabel('Distance from terminus (km)','fontsize',16); ylabel('Elevation (m)','fontsize',16);
        subplot(subV_yr);
        set(subV_yr,'xlim',xlims,'xtick',xticks,'xticklabels',[]);
        pos = get(subV_yr,'position'); set(subV_yr,'position',[pos(1)+0.02 pos(2)-0.05 pos(3) pos(4)+0.05]); %slightly stretch y-axis
        ylabel('Speed (m/d)','fontsize',16);
        subplot(subV_mo);
        set(subV_mo,'xlim',xlims,'xtick',xticks,'xticklabels',xticks/1000);
        pos = get(subV_mo,'position'); set(subV_mo,'position',[pos(1)+0.02 pos(2) pos(3) pos(4)+0.05]); %slightly stretch y-axis
        xlabel('Distance from terminus (km)','fontsize',16); ylabel('Speed (m/d)','fontsize',16);
        subplot(subT);
        pos = get(subT,'position'); set(subT,'position',[pos(1)+0.02 pos(2) pos(3) pos(4)]);
        drawnow;
        saveas(gcf,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-centerline-elev-speed-terminus_subplots.png'],'png'); %save the plots

        %clear profile variables
        clear im im_subset LCdir zprofs zdate* mean_prof seaward_idx centerline C vel_pts term_trace melmask* yrs mos zyrs zmos;
        clear DEM_num pm pt py sub* tran_* *xlims xticks *ylims yticks term_ref *_dist *_reldist *pos;
        clear bergdist dts inland_idx seaward_idx Tdate* term_* Tmos Tyrs;

        %save the structure with the centerline data
        save([root_dir,'GrIS-melange-characteristics.mat'],'MP','-v7.3');

        %create "climatologies" of seasonal size distributions for each site
        D = readtable([root_dir,sitenames(j,:),'/',sitenames(j,:),'-iceberg-distribution-timeseries.csv']);
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
            ps(p) = loglog(NaN,NaN,'-','color',seas_cmap(p,:),'linewidth',1.5); hold on;
            mo_ref = find(berg_mo==MP(j).D.months(p,1) | berg_mo==MP(j).D.months(p,2) | berg_mo==MP(j).D.months(p,3));
            if ~isempty(mo_ref)
                %plot the individual distributions colored by season
                subplot(subd2);
                for l = mo_ref
                    %plot the full distribution profile for each date
                    loglog(MP(j).D.area,berg_normdist(l,:),'-','color',seas_cmap(p,:),'linewidth',1.5); hold on;

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
                loglog(MP(j).D.area,MP(j).D.bergs(p,:),'-','color',seas_cmap(p,:),'linewidth',2); hold on;

            end
            clear mo_ref;
        end
        seas_leg = legend(ps,'DJF','MAM','JJA','SON');
        set(subd1,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',16);
        xlabel('Surface area (m^2)','fontsize',16); ylabel('Normalized iceberg count','fontsize',16);
        subplot(subd1); grid on; drawnow;
        set(subd2,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',16);
        xlabel('Surface area (m^2)','fontsize',16);
        subplot(subd2); grid on; drawnow;
        saveas(dist_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-seasonal-iceberg-distribution_loglog.png'],'png'); %save the plot
        clear D berg_* ps ;

        %save the structure with the seasonal distribution data
        save([root_dir,'GrIS-melange-characteristics.mat'],'MP','-v7.3');
        close all;

    end

end



%% create terminus position plots
for j = 1:length(MP)
    site_abbrev = MP(j).name; disp(site_abbrev);
    close all; drawnow;
    ts_fig = figure; set(ts_fig,'position',[50 50 600 1200]);

    %extract the transect_inc from the files (only works if only one
    %version is saved for each site but it's needed to account for adaptive
    %transect spacing depending on site size)
    shp_files = dir([root_dir,site_abbrev,'/shapefiles/',site_abbrev,'*.shp']);
    for k = 1:length(shp_files)
        if contains(shp_files(k).name,['transects_'])
            transect_inc = str2num(shp_files(k).name(end-8:end-5));
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
        try 
            Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        catch
            MP(j).T.date{p} = string(str2num(MP(j).T.date{p})-1);
            Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        end
        Tdate(p) = convert_to_decimaldate(char(MP(j).T.date{p}));
        Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        Tyrs(p) = year(Tdatetime(p)); Tmos(p) = month(Tdatetime(p));
        dists = sqrt((MP(j).T.termX(1,p)-C.X).^2 + (MP(j).T.termY(1,p)-C.Y).^2);
        [~,term_ind] = min(dists);
        % MP(j).T.termdist(1,p) = centerline_dist(term_ind);

        clear dists term_ind;
    end
    MP(j).T.termdist(Tyrs < min(years) | Tyrs > max(years)) = NaN;

    %compile the data for identifying terminus gaps
    term_decidates = [Tdate, zdate(term_trace ==1)]; term_dates = [Tdatetime, zdatetime(term_trace ==1)]; 
    term_X = [MP(j).T.termX, MP(j).Z.termX(term_trace ==1)];
    term_Y = [MP(j).T.termY, MP(j).Z.termY(term_trace ==1)];
    term_dists = [MP(j).T.termdist, MP(j).Z.termdist(term_trace ==1)];

    %plot the image-based termini first, then DEM
    for p = 1:length(MP(j).T.date)
        %add dummy plot for the legend
        if p ==1
            for k = 1:4
            pterm(k) = plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).T.termdist(p),Tdate(p),'x',...
            'color',seas_cmap(k,:),'linewidth',2,'markersize',10); hold on;
            end
        end

        %plot the data
        plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).T.termdist(p),Tdate(p),'x',...
            'color',seas_cmap(floor(Tmos(p)/4)+1,:),'linewidth',2,'markersize',10); hold on;
    end
    for p = 1:length(MP(j).Z.date)
        if term_trace(p) == 0 && isnan(MP(j).Z.termdist(p))
            % MP(j).Z.termX(p) = NaN; MP(j).Z.termY(p) = NaN; MP(j).Z.termdist(p) = NaN;

            % %fill gaps with data
            % dts = 365*(zdate(p)-term_decidates);
            % [mindt,minind] = min(abs(dts));
            % if mindt < 30
            %     MP(j).Z.termX(p) = term_X(minind); MP(j).Z.termY(p) = term_Y(minind);
            %     MP(j).Z.termdist(p) = term_dists(minind);
            %     plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).Z.termdist(p),zdate(p),'s',...
            %         'color','k','linewidth',1,'markerfacecolor','none','markersize',10); hold on;
            %     disp(['Filled terminus for ',MP(j).name,' ',char(MP(j).Z.date(p)),' with ',char(term_dates(minind))]);
            % else
                disp(['Need terminus data for ',MP(j).name,' ',char(MP(j).Z.date(p))]);
            % end
            % clear mindt minind;
        else
             plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).Z.termdist(p),zdate(p),'s',...
                 'color',seas_cmap(floor(zmos(p)/4)+1,:),'linewidth',1,'markerfacecolor','none','markersize',10); hold on;
        end
       
    end
    set(gca,'xlim',[0,max([MP(j).T.termdist,MP(j).Z.termdist])-min([MP(j).T.termdist,MP(j).Z.termdist])],'ylim',[min(years) max(years)]);
    xlims = get(gca,'xlim'); xticks = get(gca,'xtick');
    set(gca,'xtick',xticks,'xticklabels',xticks/1000,'fontsize',16);
    term_leg = legend(pterm,season_names,'Location','northoutside',...
            'Orientation','horizontal');
    xlabel('Distance from most-retreated terminus (km)','fontsize',16); ylabel('Year','fontsize',16);
    grid on; drawnow;
    
    %save the figure & updated terminus positions
    % save([root_dir,'GrIS-melange-characteristics.mat'],'MP','-v7.3');
    saveas(ts_fig,[root_dir,MP(j).name,'/',MP(j).name,'-terminus-timeseries-plot.png'],'png'); %save the plots
    % uiwait %advance only after figure is closed

    clear C centerline_dist dts inland_idx seaward_idx mean_prof Tdate* term_* Tmos tran_dist Tyrs zdate* zmos zyrs;
end
close all;
disp('Done plotting terminus timeseries');


%% extract melange attributes & estimate buttressing
close all; drawnow;

%REDO PLOTS TO HAVE 3 COLUMNS & 7 ROWS(?), LEAVING A BIG HOLE IN THE MIDDLE
%FOR THE MAP AND A GAP ALONG THE RIGHT EDGE FOR THE EAST COAST
% rows = 9; cols = 2; %9 sites in west Greenland
% plot_locs = [1,3,5,7,9,11,13,15,17,18,16,14,12,2];
rows = 7; cols = 3; %9 sites in west Greenland
plot_locs = [2,1,4,7,10,13,16,19,20,21,18,15,12,3];

%Thickness parameters:
zcutoff = 3; %elevation threshold below which to ignore iceberge (m)
rho_i = 900; %ice density (kg/m^3)
rho_w = 1026; %water density (kg/m^3)
Hcutoff = round((rho_w/(rho_w-rho_i))*zcutoff); %H threshold for figure naming

%Size distribution parameters
nthresh = 1e-6; % set small number bin cutoff (n1 must be greater than this value)
zthresh = 3; %set small size bin cutoff (freeboard must exceed this value)
rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)
ARcomp.best.autoALL = 2; % iceberg aspect ratio 
vthresh = (1/4)*pi*((rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL.*zthresh).^2; %don't include this bin size or smaller in curve fitting
dplawthresh = 10^5; % upper bound on the intercept for the dummy powerlaw
norm_type = 2; % toggle between L2, max, and log norm using 2, Inf, and 'log'
normalize_exp = 1.5; % Increase to weight residuals towards end of the curve, minimum = 1

%Velocity parameters: ideally only use velocities from time periods with
%DEMs to create the seasonal profiles but you can also use all velocities
%from that season (even if there is no DEM for that season-year combination)
vfilter = 'annual'; %options: 'annual' OR 'all'
vdtmin = 15; %maximum image acquisition separation for velocities
vdtmax = 60; %maximum image acquisition separation for velocities
dvterm = 'yes'; %include velocities from the terminus

%decide whether to sample data at (1) fixed locations or (2) relative
%locations with respect to the terminus positions in the DEMs
sampling = 'dated'; % options: 'dated' OR 'fixed';

% create summary figures
Hfig = figure; set(Hfig,'position',[50 50 1200 1200]);
% Vfig = figure; set(Vfig,'position',[150 50 1200 1200]);
for j = 1:length(plot_locs)
    figure(Hfig);
    eval(['subH',num2str(geo_ind(j)),'=subplot(',num2str(rows),',',num2str(cols),',',num2str(plot_locs(j)),');']);
    % figure(Vfig);
    % eval(['subV',num2str(geo_ind(j)),'=subplot(',num2str(rows),',',num2str(cols),',',num2str(plot_locs(j)),');']);
end
missfig = figure; set(missfig,'position',[950 50 600 500]); 
% subm1 = subplot(2,1,1); subm2 = subplot(2,1,2);

% %If adding new dates, profiles need to be removed from the structure
% %because the inland reference may have changed
% %(will throw errors when concatenating)
% for j = 1:length(MP)
%     MP(j).Z = rmfield(MP(j).Z,{'dist','Hseas'});
%     MP(j).V = rmfield(MP(j).V,{'dist','Vseas','dVdx'});
% end

disp(['Creating profile plots ignoring icebergs thinner than ',num2str(Hcutoff),'m & using ',vfilter,' speeds w/ dt ~',num2str(round([vdtmin,vdtmax])),' days'])

%iterate
for j = 1:length(MP)
    sitenames(j,:) = MP(j).name; seasons = MP(j).D.months;
    disp(sitenames(j,:)); site_abbrev = MP(j).name;
    transect_inc = round(nanmean(sqrt((diff(MP(j).V.X).^2 + diff(MP(j).V.Y).^2)))/100)*100;
    
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
        % seaward_idx = inland_idx-(inland_meanidx-seaward_meanidx);
        seaward_idx = inland_idx-floor(nanmean(inland_idx-seaward_ext));
    else
        error('Profile sampling strategy is not properly defined: sampling must be ''fixed'' or ''dated''')
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
                    z_temp = flipud(Zfilt(1:inland_idx(yr_idx(k))-1,yr_idx(k)));
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
                pack_seas(:,k,p) = nanmean(pack_profiles(ismember(mos_yr,seasons(k,:))==1,:),1)';
                if ~isempty(inland_idx(ismember(mos_yr,seasons(k,:))==1))
                    MP(j).B.ref(1,k,p) = mean(inland_idx(ismember(mos_yr,seasons(k,:))==1))-1;
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
    z_seas(z_seas==0) = NaN; H_seas(H_seas==0) = NaN; pack_seas(pack_seas==0) = NaN; 
    zdist = (transect_inc/2):transect_inc:(max(inland_idx)-1)*transect_inc; 
    Hdist = zdist-transect_inc/2; MP(j).Z.dist = Hdist;

    %if using all velocities for seasonal profiles, not just from DEM
    %years, create a timeseries of the inland_idx based off the terminus
    %delineations
    for k = 1:length(MP(j).T.termX)
        try 
            Tdatetime(k) = datetime(MP(j).T.date{k},'InputFormat','yyyyMMdd');
        catch
            MP(j).T.date(k) = string(str2num(MP(j).T.date{k})-1);
            Tdatetime(k) = datetime(MP(j).T.date{k},'InputFormat','yyyyMMdd');
        end
        Tdate(k) = convert_to_decimaldate(char(MP(j).T.date{k}));
        Tdatetime(k) = datetime(MP(j).T.date{k},'InputFormat','yyyyMMdd');
        Tyrs(k) = year(Tdatetime(k)); Tmos(k) = month(Tdatetime(k));
        dists = sqrt((MP(j).T.termX(k)-MP(j).V.X).^2 + (MP(j).T.termY(k)-MP(j).V.Y).^2);
        [~,ia] = sort(dists);
        term_idx(k) = max(ia(1:2));
        clear dists ia;
    end
    for p = 1:length(years)
        yr_idx = find(Tyrs == years(p));
        vterm_idx(p) = round(median(term_idx(yr_idx)));
    end

    %grab velocity timeseries for the transect-centerline intersection
    %points & create seasonally-averaged climatologies
    if contains(dvterm,'y')
        inland_vel = max(inland_idx); %grab data from the terminus
        ref_adjust = 1; %shift the referencing for the relative coordinate system
    else
        inland_vel = max(inland_idx)-1; %only get speeds from the melange
        ref_adjust = 0;
    end
    %loop through the pts
    vel_pts = dir([root_dir,sitenames(j,:),'/velocities/']);
    vel_seas = NaN(max(inland_idx),4,length(years));
    for i = 1:length(vel_pts)
        if contains(vel_pts(i).name,'velocity')
            pt_ref = str2num(vel_pts(i).name(end-5:end-4));

            if pt_ref <= inland_vel %setting this <= and making rel_ref = inland_idx - pt_ref+1 gives the inland pt on the glacier

                %read the file
                V = readtable([root_dir,sitenames(j,:),'/velocities/',vel_pts(i).name]);

                %filter out all the velocities based on temporal resolution
                short_dts = find(V.days_dt>vdtmin & V.days_dt<vdtmax); %get rid of all velocities with coarse temporal resolution
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
                                    rel_ref = vterm_idx(p)-pt_ref+ref_adjust; %edit if monthly terminus positions are available
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
                                        rel_ref = round(median([inland_idx(zyrs==years(p) & ismember(zmos,seasons(k,:))), term_idx(Tyrs==years(p) & ismember(Tmos,seasons(k,:)))]))-pt_ref+ref_adjust; 
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
    vdist = (transect_inc*[0:1:size(vel_seas,1)-1])+((transect_inc/2)-(ref_adjust*transect_inc)); %-(transect_inc/2) start for inland pt on the glacier, (transect_inc/2) start for inland point in the melange
    MP(j).V.dist = vdist;

    %loop through the subsetted size distributions and compute seasonal
    %averages nearest the terminus & a fixed relative distance down-fjord
    Dsubs = dir([root_dir,MP(j).name,'/',MP(j).name,'*-iceberg-distribution-subsets.csv']);
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(j).name,'/',Dsubs(p).name],"VariableNamingRule","preserve");
        berg_datestrings(p,:) = Dsubs(p).name(length(site_abbrev)+2:length(site_abbrev)+9);
        % yrs(p) = year(berg_dates(p));
        mos(p) = month(berg_dates(p));

        if ~isnan(inland_idx(p))
            %full size distributions at two points
            bergdist_inland(p,:) = table2array(D(:,(inland_idx(p)-1)+2))'; %near-terminus
            bergdist_seaward_setdx(p,:) = table2array(D(:,seaward_meanidx+2))'; %standard distance from near-terminus
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
    %load or create iceberg distribution profiles
    if exist([root_dir,MP(j).name,'/models/',MP(j).name,'-powerlaw-slope_seasonal-profiles.csv']) == 2
        D = readtable([root_dir,MP(j).name,'/models/',MP(j).name,'-powerlaw-slope_seasonal-profiles.csv'],"VariableNamingRule","preserve");
        size_plslope = table2array(D(:,2:end)); clear D;
        D = readtable([root_dir,MP(j).name,'/models/',MP(j).name,'-powerlaw-bergybit-mispredict_seasonal-profiles.csv'],"VariableNamingRule","preserve");
        res_A = table2array(D(:,2:end)); clear D;
        %find the last binned seasonal distribution with data
        for k = 1:4
            if sum(~isnan(res_A(k,:))) > 0
                res_Alim(k) = find(~isnan(res_A(k,:))==1,1,'last');
            else
                res_Alim(k) = NaN;
            end
        end
    else
        disp('Fitting empirical curves to binned, seasonal size distributions:')
        Dsubs = dir([root_dir,MP(j).name,'/',MP(j).name,'*-iceberg-distribution-subsets.csv']);
        for p = 1:length(Dsubs)
            D = readtable([root_dir,MP(j).name,'/',Dsubs(p).name],"VariableNamingRule","preserve");

            %initialize the data compilation
            if p == 1
                v1 = table2array(D(:,1)); dv1 = table2array(D(:,2));
                n1_seas = NaN(length(v1),max(inland_idx),4);
            end

            %first 2 columns are area & bin width, so 3+ are data from points along the centerline
            n1 = NaN(length(v1),max(inland_idx)-1);
            if ~isnan(inland_idx(p))
                berg_nos = fliplr(table2array(D(:,3:(inland_idx(p)-1)+2)));
                berg_nos(berg_nos==0) = NaN;
                n1(:,1:size(berg_nos,2)) = berg_nos./dv1;
            end
            clear D;

            %stack distributions in a structure
            binned_n1_stack(:,:,p) = n1;
            clear n1 berg_nos;
        end
        %create seasonal average distributions for each fjord bin
        for k = 1:4
            disp([char(season_names{k}),': parameters starting from the terminus & moving seaward...']);
            
            %stack the profiles
            n1_seasmean(:,:,k) = nanmean(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),3);
            n1_seasmax(:,:,k) = (nanmean(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),3)+std(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),0,3,'omitnan'));
            n1_seasmin(:,:,k) = (nanmean(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),3)-std(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),0,3,'omitnan'));

            for l = 1:size(n1_seasmean,2)
                n1 = n1_seasmean(:,l,k);

                %filter to decide if there are enough points to fit lines
                v = v1(n1>nthresh); n = n1(n1>nthresh); dv = dv1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
                n = n(v>vthresh); dv = dv(v>vthresh); v = v(v>vthresh); %remove smallest size classes
                v = v(~isnan(n)); dv = dv(~isnan(n)); n = n(~isnan(n)); % Remove NaNs

                % Fit E-BC fragmentation model to the data and plot
                if sum(~isnan(v)) > 12
                    %enter unfiltered data into model (filtering happens there)
                    [alpha,c1,c2,c3,c4,~,~] = EBC_fragmentation_curve(v1, n1, dv1, norm_type, vthresh, nthresh, 1, normalize_exp); % fit Eq. 1
                    n1_mod = EBC_model([c1,c2,alpha,c3,c4],v1); % grab model

                    %calculate residuals: negative residuals indicate that
                    %the powerlaw fit under-estimates icebergs
                    res = n1_mod-n1; % calculate residuals from E-BC model (data - model)
                    % resp = res(res >= 0); resn = res(res < 0); % split pos. and neg. residuals
                    % resp_idx = find(res >= 0); % grab index of positive residuals
                    % resn_idx = find(res < 0); % grab index of negative residuals

                    % save the parameters
                    size_plslope(k,l) = alpha; res_no(k,l) = sum(res(v1<=vthresh).*dv1(v1<=vthresh)); 
                    res_A(k,l) = sum(v1(v1<=vthresh).*(res(v1<=vthresh).*dv1(v1<=vthresh)))./sum(v1.*(n1.*dv1),"omitnan");
                    % clear variables
                    clear res n_mod v dv n n1 alpha c1 c2 c3 c4;
                else
                    size_plslope(k,l) = NaN; res_no(k,l) = NaN; res_A(k,l) = NaN;
                    clear v dv n n1;
                end
            end

            %find the last binned seasonal distribution with data
            if sum(~isnan(res_A(k,:))) > 0
                res_Alim(k) = find(~isnan(res_A(k,:))==1,1,'last');
            else
                res_Alim(k) = NaN;
            end
        end
        Tslopes = [table(season_names'), array2table(size_plslope)]; Tres = [table(season_names'), array2table(res_A)];
        for l = 1:size(size_plslope,2)
            bin_dist(l) = 0+(l-1)*transect_inc;
            bin_slopename(l) = {['Powerlaw exp: ',num2str(bin_dist(l)),'m']};
            bin_residname(l) = {['Bergybit misprediction: ',num2str(bin_dist(l)),'m']};
        end
        column_names = ["Season", bin_slopename]; Tslopes.Properties.VariableNames = column_names; 
        column_names = ["Season", bin_residname]; Tres.Properties.VariableNames = column_names;
        writetable(Tslopes,[root_dir,MP(j).name,'/models/',MP(j).name,'-powerlaw-slope_seasonal-profiles.csv']);
        writetable(Tres,[root_dir,MP(j).name,'/models/',MP(j).name,'-powerlaw-bergybit-mispredict_seasonal-profiles.csv']);
    end
    disp('Done extracting size distribution information')

    %Site-specific plots
    sitefig = figure; set(sitefig,'position',[50 850 1200 1000]);
    subz = subplot(4,2,[1:2]); subv = subplot(4,2,[3:4]); 
    max_xlim = max([Hdist(find(sum(sum(~isnan(H_seas),2),3)>0,1,'last'))]);
    plot_ind = find(geo_ind == j);

    %add thickness profiles to figures
    for k = 1:4
        if k <= 2
            fill_alpha = 0.1;
        else
            fill_alpha = 0.2;
        end

        %relative to moving terminus
        Hmean = nanmean(H_seas(:,k,:),3);
        Hmax = (nanmean(H_seas(:,k,:),3)+std(H_seas(:,k,:),0,3,'omitnan')); 
        Hmin = (nanmean(H_seas(:,k,:),3)-std(H_seas(:,k,:),0,3,'omitnan')); 

        %plot
        Hmax_idx = find(~isnan(Hmax)==1); Hmin_idx = find(~isnan(Hmin)==1);
        if sum(~isnan(Hmean)) > 0
            figure(sitefig); subplot(subz);
            % fill([zdist(Hmax_idx), fliplr(zdist(Hmin_idx))]',[Hmax(Hmax_idx); flipud(Hmin(Hmin_idx))],...
            %     seas_cmap(k,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
            % pz(k) = plot(zdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
            fill([Hdist(Hmax_idx), fliplr(Hdist(Hmin_idx))]',[Hmax(Hmax_idx); flipud(Hmin(Hmin_idx))],...
                seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
            pz(k) = plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;

            %regional fig
            figure(Hfig);
            if ~isempty(plot_ind)
                %navigate to the subplot
                eval(['subplot(subH',num2str(j),');']);
                ax1 = gca;

                %plot the thickness profiles
                yyaxis left
                fill([Hdist(Hmax_idx), fliplr(Hdist(Hmin_idx))]',[Hmax(Hmax_idx); flipud(Hmin(Hmin_idx))],...
                    seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
                pH(k) = plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
            end
            MP(j).Z.Hseas(k,:) = Hmean';
        else
           figure(sitefig); subplot(subz);
           pz(k) = plot(NaN,NaN,'-','color',seas_cmap(k,:),'linewidth',3); hold on; 
           MP(j).Z.Hseas(k,:) = NaN(1,size(H_seas,1));
        end
        % H_ylim(k,:) = [min(Hmin), max(Hmax)];
        H_ylim(k,:) = [min(Hmean), max(Hmean)];
        clear Hmean Hmax* Hmin*;
        % clear zdist;
    end
    %format the thickness subplot in the site figure
    figure(sitefig); subplot(subz);
    seas_leg = legend(pz,season_names);
    set(gca,'fontsize',16); drawnow; grid on; 
    set(subz,'xlim',[0,max_xlim],...
        'xticklabel',[],'ylim',[floor(min(H_ylim(:,1))/10)*10 ceil(max(H_ylim(:,2))/10)*10]); 
    ylims = get(gca,'ylim'); 
    if range(ylims) <= 50
        set(gca,'ytick',[floor(min(H_ylim(:,1))/10)*10:10:ceil(max(H_ylim(:,2))/10)*10]);
    else
        set(gca,'ytick',[floor(min(H_ylim(:,1))/10)*10:20:ceil(max(H_ylim(:,2))/10)*10]);
    end
    % set(subz,'xlim',[0,floor(tran_reldist(seaward_idx)/1000)*1000+500],'xticklabel',[],'ylim',[0 ceil(max(H_ylim)/50)*50]); 
    clear ylims;
    ylabel('Thickness (m)','fontsize',16); %xlabel('Distance from terminus (km)','fontsize',16); 
    subz_pos = get(subz,'position'); 
    seas_leg.Location = 'northoutside'; seas_leg.Orientation = 'horizontal'; %move the legend
    set(subz,'position',subz_pos); %resize the subplot
    %format the site subplot in the thickness figure
    if ~isempty(plot_ind)
        figure(Hfig); eval(['subplot(subH',num2str(j),');']);
        yyaxis left;
        set(gca,'xlim',[0, 25000],'xticklabel',[],...
            'ylim',[0 ceil(max(H_ylim(:,2))/10)*10],'fontsize',16); drawnow;
        gca_pos = get(gca,'position'); ylims = get(gca,'ylim'); 
        if max(ylims) <= 55
            set(gca,'ylim',[-10,55],'ytick',[25:10:55]);
        else
            set(gca,'ylim',[-50,100],'ytick',[25:25:100]);
        end
        % yticks = get(gca,'ytick'); set(gca,'ytick',yticks(2:end),'yticklabel',yticks(2:end));
        ax1.YAxis(1).Color = 'k'; plot([0, 25000],[25,25],'-k','linewidth',1.5);
        set(ax1,'box','on'); ax1.LineWidth = 1.5;
        %add labels
        if plot_locs(plot_ind) >= rows*cols - (cols-1) || (plot_locs(plot_ind) >= 2 && plot_locs(plot_ind) <= 3)
            xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks/1000); clear xticks;
            % set(gca,'yticklabel',yticks);
            xlabel('Distance from terminus (km)','fontsize',16); 
            if plot_locs(plot_ind) == rows*cols - (cols-1)
                ax1.YAxis(1).Label.String = 'Thickness (m)';
                ax1.YAxis(1).Label.Position = [-2.7066e3 260 -1];
                % ylabel('Thickness (m)','fontsize',16);
            end
        end
        clear ylims yticks;
        %add a legend
        if j == length(MP)
            H_leg = legend(pH,season_names);
            subH_pos = get(gca,'position');
            H_leg.Orientation = 'horizontal'; %move the legend
            H_leg.Position = [0.39 0.24 0.25 0.035]; 
            set(gca,'position',subH_pos); %resize the subplot
            drawnow;
        end
        set(gca,'position',[gca_pos(1) gca_pos(2) gca_pos(3) 0.07]);
    end
    
    %add speed profiles to figures
    for k = 1:4
        if k <= 2
            fill_alpha = 0.1;
        else
            fill_alpha = 0.2;
        end

        %relative to moving terminus
        vmean = nanmean(vel_seas(:,k,:),3)./365;
        vmax = (nanmean(vel_seas(:,k,:),3)+std(vel_seas(:,k,:),0,3,'omitnan'))./365; 
        vmin = (nanmean(vel_seas(:,k,:),3)-std(vel_seas(:,k,:),0,3,'omitnan'))./365; 
        
        %mask out values near zero because they are so much slower than the
        %melange velocities that they must be bad sea ice tracking
        bottomout = find(vmean < 0.5,1,'first');
        vmean(bottomout:end) = NaN; vmax(bottomout:end) = NaN; vmin(bottomout:end) = NaN; 

        %plot
        vmax_idx = find(~isnan(vmax)==1); vmin_idx = find(~isnan(vmin)==1);
        if sum(~isnan(vmean)) ~= 0
            figure(sitefig); subplot(subv);
            fill([vdist(vmax_idx), fliplr(vdist(vmin_idx))]',[vmax(vmax_idx); flipud(vmin(vmin_idx))],...
                seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
            pv(k) = plot(vdist(~isnan(vmean))',vmean(~isnan(vmean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;

            %regional fig
            figure(Hfig);
            % figure(Vfig);
            if ~isempty(plot_ind)
                %navigate to the subplot
                eval(['subplot(subH',num2str(j),');']);
                % eval(['subplot(subV',num2str(j),');']);

                %plot the thickness profiles
                yyaxis right;
                fill([vdist(vmax_idx), fliplr(vdist(vmin_idx))]',[vmax(vmax_idx); flipud(vmin(vmin_idx))],...
                    seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
                if k >= 2
                    pV(k-1) = plot(vdist(~isnan(vmean))',vmean(~isnan(vmean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
                else
                    plot(vdist(~isnan(vmean))',vmean(~isnan(vmean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
                end
            end
            MP(j).V.Vseas(k,:) = vmean';
        else
            figure(sitefig); subplot(subv);
            MP(j).V.Vseas(k,:) = NaN(1,size(vel_seas,1));
        end
        clear vmean vmax* vmin*;
        % clear vdist;
    end
    %format the speed subplot in the site figure
    figure(sitefig); subplot(subv);
    set(gca,'fontsize',16); grid on; drawnow;
    set(subv,'xlim',[0,max_xlim]);
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
    xlabel('Distance from terminus (km)','fontsize',16); ylabel('Speed (m/d)','fontsize',16);
    pos = get(subv,'position'); %set(subv,'position',[pos(1) pos(2)+0.05 pos(3) pos(4)]);
    %format the site subplot in the speed figure
    if ~isempty(plot_ind)
        figure(Hfig); eval(['subplot(subH',num2str(j),');']); yyaxis right;
        % figure(Vfig); eval(['subplot(subV',num2str(j),');']);
        set(gca,'xlim',[0,25000],'xticklabel',[],'fontsize',16); drawnow;
        % gca_pos = get(gca,'position');
        ylims = get(gca,'ylim'); set(gca,'ylim',[0 max(ylims)]);
        if max(ylims) <= 35
            set(gca,'ylim',[0,70],'ytick',[0:10:30]);
        % elseif max(ylims) <= 40
        %     set(gca,'ylim',[0,40],'ytick',[0:20:40]);
        else
            set(gca,'ylim',[0,140],'ytick',[0:20:60]);
        end
        yticks = get(gca,'ytick'); 
        ax1.YAxis(2).Color = 'k'; 
        if plot_locs(plot_ind) >= rows*cols - (cols-1)
            xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks/1000); clear xticks;
            set(gca,'yticklabel',yticks);
            xlabel('Distance from terminus (km)','fontsize',16);
            if plot_locs(plot_ind) == rows*cols
                ax1.YAxis(2).Label.String = 'Speed (m/d)';
                ax1.YAxis(2).Label.Position = [2.7559e4 530 -1];
                % ylabel('Speed (m/d)','fontsize',16);
            end
        end
        clear ylims yticks;
        % %add a legend
        % figure(Vfig); eval(['subplot(subV',num2str(j),');']);
        % if j == length(MP)
        %     V_leg = legend(pV,season_names(2:end));
        %     subV_pos = get(gca,'position');
        %     V_leg.Orientation = 'horizontal'; %move the legend
        %     V_leg.Position = [0.4 0.945 0.2 0.035]; 
        %     set(gca,'position',subV_pos); %resize the subplot
        %     drawnow;
        % end
        % set(gca,'position',[gca_pos(1) gca_pos(2) gca_pos(3) 0.07]);
    end
    gca_pos = get(gca,'position');
    set(gca,'position',[gca_pos(1) gca_pos(2) gca_pos(3) 1.5*gca_pos(4)]);
    drawnow; clear plot_ind pv pV;
    
    %add subplots beneath the figure that generally show the location of the size distributions
    %inland bin
    figure(sitefig); subplot(4,2,5);
    for k = 1:4
        loglog(MP(j).D.area,bergdist_seas(k,:,2),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
    end
    set(gca,'ylim',[10^-10 10^-1],'xlim',[10^1 10^7],...
        'ytick',[10^-9,10^-5,10^-1],...
        'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',16); grid on;
    xlabel('Surface area (m^2)','fontsize',16); ylabel('Iceberg count','fontsize',16);
    text(10000,10^-2,'near-terminus','fontsize',16);
    %seaward bin
    subplot(4,2,6);
    for k = 1:4
        loglog(MP(j).D.area,bergdist_seas(k,:,1),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
    end
    set(gca,'ylim',[10^-10 10^-1],'xlim',[10^1 10^7],...
        'ytick',[10^-9,10^-5,10^-1],...
        'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',16); grid on;
    xlabel('Surface area (m^2)','fontsize',16);
    text(10000,10^-2,'seaward','fontsize',16);
    drawnow;
    
    %estimate buttressing
    if contains(vfilter,'annual')
        %seasonal averages
        for k = 1:4
            %add seasonal strain rate profiles to the structure
            dvdist_temp = diff(vdist);
            MP(j).V.dVdx(k,:) = (diff(MP(j).V.Vseas(k,:))./dvdist_temp)/365;
            sr_pivot = find(MP(j).V.dVdx(k,:)<0,1,'first');
            if ~isempty(sr_pivot)
                sr_dsign_pt(k) = vdist(sr_pivot+ref_adjust);
            else
                sr_dsign_pt(k) = NaN;
            end
            clear sr_pivot;

            %calculate buttressing for each year
            for p = 1:length(years)
                % add seasonal buttressing info to the structure
                MP(j).B.dVdx(1,k,p) = (diff(vel_seas(1+ref_adjust:2+ref_adjust,k,p))./diff(vdist(1+ref_adjust:2+ref_adjust)))/365; %shift referencing from velocities as needed
                press = 0.5*rho_i*(1-(rho_i/rho_w))*9.81*MP(j).B.Ho(zcutoff+1,k,p);
                MP(j).B.butt_Meng(zcutoff+1,k,p) = press*MP(j).B.packing(zcutoff+1,k,p)*MP(j).B.Ho(zcutoff+1,k,p);
                MP(j).B.butt_Amundson(zcutoff+1,k,p) = (-2*(MP(j).B.Ho(zcutoff+1,k,p)*press*MP(j).B.dVdx(1,k,p))/((MP(j).B.dVdx(1,k,p)/0.3)+MP(j).B.dVdx(1,k,p)))+press*MP(j).B.Ho(zcutoff+1,k,p);
                clear press dvdist_temp;
            end
        end

        %CALCULATE BUTTRESSING FOR EACH INDIVIDUAL DATE 
        % (to create probability density functions showing
        %the full suite of buttressing estimates)


    end

    %plot seasonal profiles of melange characteristics
    % bergAfig = figure; set(gcf,'position',[50 50 1200 400]); 
    figure(sitefig); subplot(4,2,[7,8]);
    clear pz;
    for k = 1:4
        %site figure
        % figure(bergAfig);
        figure(sitefig); 
        ax2 = gca;
        yyaxis left
        % add dummy data for legend beneath real data
        if sum(~isnan(size_plslope(k,:))) > 0
            for l = 1:4
                ps(l) = plot(Hdist,res_A(k,:),'-','color',seas_cmap(l,:),'linewidth',3); hold on; %powerlaw fits
            end
        end

        %plot the fractional difference between powerlaw-modeled & observed
        %area covered by small icebergs
        plot(Hdist,res_A(k,:),'-','color',seas_cmap(k,:),'linewidth',3); hold on; %bergy bit misprediction
        ax2.YLim = [-1,1];

        %plot packing density
        Hmean = nanmean(H_seas(:,k,:),3);
        packmean = 100*nanmean(pack_seas(:,k,:),3);
        if sum(~isnan(Hmean)) > 0
            yyaxis right;
            % pz(k) = plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
            % yyaxis right;
            plot(Hdist(~isnan(Hmean))',packmean(~isnan(Hmean)),'--','color',seas_cmap(k,:),'linewidth',3); hold on;
            ax2.YLim = [0,100]; 
        end

        %format the site figure
        ax2.YAxis(1).Color = 'k'; ax2.YAxis(2).Color = 'k'; 
        ax2.YAxis(1).Label.String = 'Bergy bit misprediction';
        ax2.YAxis(2).Label.String = 'Packing density (%)';

        %scatterplot of bergy bit misfits vs packing density for all sites
        figure(missfig); set(gca,'box','on');
        if ismember(MP(j).name,big3)
            %plot all data along the profiles
            scatter(res_A(k,:),packmean,200-sqrt(Hdist),'d','MarkerFaceColor','none',...
                'MarkerEdgeColor',seas_cmap(k,:),'LineWidth',1.5); hold on;

            % %plot only data from at the terminus & seaward margins
            % %inland
            % subplot(subm1);
            % plot(res_A(k,1),packmean(1),'d','MarkerFaceColor','none',...
            %     'MarkerEdgeColor',seas_cmap(k,:),'MarkerSize',13,'LineWidth',1.5); hold on;
            % 
            % %seaward
            % subplot(subm2);
            % plot(res_A(k,max(inland_idx)-seaward_meanidx+1),packmean(max(inland_idx)-seaward_meanidx+1),'d','MarkerFaceColor','none',...
            %     'MarkerEdgeColor',seas_cmap(k,:),'MarkerSize',13,'LineWidth',1.5); hold on;
        else
            %plot all data along the profiles
            scatter(res_A(k,:),packmean,200-sqrt(Hdist),'s','MarkerFaceColor','none',...
                'MarkerEdgeColor',seas_cmap(k,:),'LineWidth',1.5); hold on;
            % plot(res_A(k,:),packmean,'s','MarkerFaceColor','none',....
            %     'MarkerEdgeColor',seas_cmap(k,:),'MarkerSize',14,'LineWidth',1.5); hold on;

            % %plot only data from at the terminus & seaward margins
            % %inland
            % subplot(subm1);
            % plot(res_A(k,1),packmean(1),'s','MarkerFaceColor','none',...
            %     'MarkerEdgeColor',seas_cmap(k,:),'MarkerSize',13,'LineWidth',1.5); hold on;
            % 
            % %seaward
            % subplot(subm2);
            % plot(res_A(k,max(inland_idx)-seaward_meanidx+1),packmean(max(inland_idx)-seaward_meanidx+1),'s','MarkerFaceColor','none',...
            %     'MarkerEdgeColor',seas_cmap(k,:),'MarkerSize',13,'LineWidth',1.5); hold on;
        end
        disp(['Min bergy bit plot symbol size: ',num2str(min(200-sqrt(Hdist)))]);
        drawnow; clear packmean Hmean;
    end
    %format the site figure
    % figure(bergAfig);
    figure(sitefig); 
    % pos = get(gca,'position');
    set(gca,'fontsize',16); set(gca,'xlim',[0,max_xlim]); 
    xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks/1000); clear xticks;
    xlabel('Distance from terminus (km)','fontsize',16); 
    yyaxis left
    plot([0,max_xlim],[0,0],'-k','linewidth',1); hold on; %add zero misfit line
    % %add a dummy legend for line types
    % rectangle(gca,'Position',[0.85*max_xlim, 0.5*max(get(gca,'ylim')), 0.14*max_xlim, 0.45*max(get(gca,'ylim'))],...
    %     'FaceColor','w','EdgeColor','k');
    % plot([0.86*max_xlim,0.89*max_xlim],[0.825*max(get(gca,'ylim')),0.825*max(get(gca,'ylim'))],'-k','linewidth',3); hold on;
    % plot([0.86*max_xlim,0.89*max_xlim],[0.625*max(get(gca,'ylim')),0.625*max(get(gca,'ylim'))],'--k','linewidth',3); hold on;
    % text(0.90*max_xlim,0.81*max(get(gca,'ylim')),'left axis','fontsize',18);
    % text(0.90*max_xlim,0.61*max(get(gca,'ylim')),'right axis','fontsize',18);
    %plot a seasonal legend
    % seas_leg = legend(ps,season_names); set(seas_leg,'location','northoutside','orientation','horizontal');
    % leg_pos = get(seas_leg,'position'); set(seas_leg,'position',[leg_pos(1) 0.94 leg_pos(3) leg_pos(4)]);
    % set(gca,'position',[pos(1) pos(2) pos(3) pos(4)+0.05]);
    %format the scatterplot
    figure(missfig);
    % subplot(subm1);
    % set(gca,'fontsize',16); grid on; drawnow;
    % set(gca,'xlim',[-1,1.5],'ylim',[0,100]);
    % ylabel('Packing density (%)','fontsize',16); 
    % subplot(subm2);
    set(gca,'fontsize',16); grid on; drawnow;
    set(gca,'xlim',[-1,1.5],'ylim',[0,100]);
    xlabel('Bergy bit misprediction (% Area)','fontsize',16); ylabel('Packing density (%)','fontsize',16); 
    clear pz ps *_leg pos;

    %display characteristics
    disp([site_abbrev,' Melange Characteristics:'])
    disp(' Near-Terminus Thickness: m')
    disp(['    winter = ',num2str(round(nanmedian(MP(j).B.Ho(zcutoff+1,1,:)),1)),', spring = ',num2str(round(nanmedian(MP(j).B.Ho(zcutoff+1,2,:)),1)),', summer = ',num2str(round(nanmedian(MP(j).B.Ho(zcutoff+1,3,:)),1)),', fall = ',num2str(round(nanmedian(MP(j).B.Ho(zcutoff+1,4,:)),1))]);
    for k = 1:4
        if ~isnan(MP(j).V.Vseas(k,1+ref_adjust))
            Vterm(k) = MP(j).V.Vseas(k,1+ref_adjust);
            dVdxterm(k) = MP(j).V.dVdx(k,1+ref_adjust);
        else
            if ref_adjust == 1
                Vterm(k) = nanmean(MP(j).V.Vseas(k,1:3));
                dVdxterm(k) = nanmean(MP(j).V.dVdx(k,1:3));
            else
                Vterm(k) = MP(j).V.Vseas(k,2);
                dVdxterm(k) = MP(j).V.dVdx(k,2);
            end
        end
    end
    disp(' Near-Terminus Speed: m/yr')
    disp(['    winter = ',num2str(round(365*Vterm(1),0)),', spring = ',num2str(round(365*Vterm(2),0)),', summer = ',num2str(round(365*Vterm(3),0)),', fall = ',num2str(round(365*Vterm(4),0))]);
    disp(' Near-Terminus Strain Rate: (m/yr)/m')
    disp(['    winter = ',num2str(round(365*dVdxterm(1),4)),', spring = ',num2str(round(365*dVdxterm(2),4)),', summer = ',num2str(round(365*dVdxterm(3),4)),', fall = ',num2str(round(365*dVdxterm(4),4))]);
    disp([' Seaward size distribution reference location = ',num2str(round(((max(inland_idx)-seaward_meanidx+1)*transect_inc)/1000,1)),'km']);
    disp(' ');
    disp(' Near-Terminus Packing Density: fractional area')
    disp(['    winter = ',num2str(round(nanmedian(MP(j).B.packing(zcutoff+1,1,:)),2)),', spring = ',num2str(round(nanmedian(MP(j).B.packing(zcutoff+1,2,:)),2)),', summer = ',num2str(round(nanmedian(MP(j).B.packing(zcutoff+1,3,:)),2)),', fall = ',num2str(round(nanmedian(MP(j).B.packing(zcutoff+1,4,:)),2))]);
    disp(' Near-Terminus Bergy Bit Overestimation: % Area')
    disp(['    winter = ',num2str(round(res_A(1,1),2)),', spring = ',num2str(round(res_A(2,1),2)),', summer = ',num2str(round(res_A(3,1),2)),', fall = ',num2str(round(res_A(4,1),2))]);
    disp(' Extension-Compression Strain Rate Change Location: km')
    disp(['    winter = ',num2str(round(sr_dsign_pt(1)/1000,2)),', spring = ',num2str(round(sr_dsign_pt(2)/1000,2)),', summer = ',num2str(round(sr_dsign_pt(3)/1000,2)),', fall = ',num2str(round(sr_dsign_pt(4)/1000,2))]);
    disp(' Change in Bergy Bit Misprediction Along Flow: % Area')
    disp(['    winter = ',num2str(round(res_A(1,min(res_Alim))-res_A(1,1),2)),', spring = ',num2str(round(res_A(2,min(res_Alim))-res_A(2,1),2)),', summer = ',num2str(round(res_A(3,min(res_Alim))-res_A(3,1),2)),', fall = ',num2str(round(res_A(4,min(res_Alim))-res_A(4,1),2))]);
    disp(' ');
    disp(' Buttressing estimated using near-terminus observations:');
    disp('   Meng Eqn (thickness- & packing-based): x10^6 N/m');
    disp(['    winter = ',num2str(round(nanmedian(MP(j).B.butt_Meng(zcutoff+1,1,:))/10^6,2)),', spring = ',num2str(round(nanmedian(MP(j).B.butt_Meng(zcutoff+1,2,:))/10^6,2)),', summer = ',num2str(round(nanmedian(MP(j).B.butt_Meng(zcutoff+1,3,:))/10^6,2)),', fall = ',num2str(round(nanmedian(MP(j).B.butt_Meng(zcutoff+1,4,:))/10^6,2))]);
    disp('   Amundson Eqn (thickness- & strainrate-based): x10^6 N/m');
    disp(['    winter = ',num2str(round(nanmedian(MP(j).B.butt_Amundson(zcutoff+1,1,:))/10^6,2)),', spring = ',num2str(round(nanmedian(MP(j).B.butt_Amundson(zcutoff+1,2,:))/10^6,2)),', summer = ',num2str(round(nanmedian(MP(j).B.butt_Amundson(zcutoff+1,3,:))/10^6,2)),', fall = ',num2str(round(nanmedian(MP(j).B.butt_Amundson(zcutoff+1,4,:))/10^6,2))]);
    disp(' ');

    %save the data and the figure
    save([root_dir,'GrIS-melange-characteristics.mat'],'MP','-v7.3');
    saveas(sitefig,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-speed-size_',num2str(Hcutoff),'m-Hthreshold_',num2str(vdtmin),'-',num2str(vdtmax),'dt-',vfilter,'-speeds_',sampling,'-profiles.png'],'png'); %save the plots
    % saveas(bergAfig,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-melange-properties_profiles.png'],'png'); %save the plots
    disp(['Done with #',num2str(j),': ',MP(j).name]);
    % figure(bergAfig);
    % disp('Close misfit/packing density figure to advance'); disp(' ');
    % uiwait %advance only after figure is closed
    
    %refresh
    clear berg_* bergdist* berg_normdist* C centerline* D Dsubs *idx seaward_* inland_* term_* tran_* size_classes Zfilt H_* Havg vel_* vels* w zdate berg_mo bins bin_no z_* pos pz pv *dist *yrs *mos seas_leg packing;
    clear leg_* ref_* Tdate* xlims Tslopes res* Tres size_plslope bin* column_names dv1 v1 n1* *_seas ax ax2 max_xlim width_prof Vterm dVdxterm subv subz*;
    clear ax1 ax2;
    close(sitefig); 
    % close(bergAfig); 
    drawnow; 
end
% %save the GrIS-wide plots
saveas(Hfig,[root_dir,'GrIS-melange-thickness_profiles.png'],'png'); 
% saveas(Vfig,[root_dir,'GrIS-melange-speed_profiles.png'],'png'); 
saveas(missfig,[root_dir,'GrIS-bergybit_scatterplot.png'],'png'); 


%% create overview maps for each site
cd(root_dir);
load([root_dir,'GrIS-melange-characteristics.mat']);

for j = 1:length(MP)
    close all; drawnow;
    disp(sitenames(j,:)); output_dir = [root_dir,sitenames(j,:),'/']; site_abbrev = MP(j).name;
    
    %load the melange masks
    cd([root_dir,sitenames(j,:)]);
    load([MP(j).name,'-melange-masks.mat']); %load the melange mask file

    %load the reference image
    ims = dir([root_dir,sitenames(j,:),'/images/',sitenames(j,:),'*S2*.tif']); 
    ref_image = [ims(1).folder,'/',ims(1).name];
    clear ims;
    %load the image geotiff
    [I,R] = readgeoraster(ref_image);
    im.x = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
    im.y = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
    if size(I,3) == 1 %panchromatic
        im.z = double(I);
    elseif size(I,3) == 3 %RGB
        im.z = double(I);
    else %all bands
        im.z = double(I(:,:,2:4));
    end
    clear I R;
    %crop the image to adjust brightnesses appropriately
    xlims = [find(im.x<=min(melmask.uncropped.x),1,'last'), find(im.x<=max(melmask.uncropped.x),1,'last')];
    ylims = [find(im.y>=max(melmask.uncropped.y),1,'last'), find(im.y<=min(melmask.uncropped.y),1,'first')];
    im_subset = im.z(min(ylims):max(ylims),min(xlims):max(xlims),:);
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

    %add manual delineations from images
    cd([root_dir,site_abbrev,'/termini/']);
    % termfiles = dir('*term*.shp');
    if contains(site_abbrev,'KBG') %skip looking for manual delineations other than termpicks for KBG
        termfile_names = {[site_abbrev,'_termpicks.shp']}; %specify consistent data sources
    else
        termfile_names = {[site_abbrev,'_termpicks.shp'],[site_abbrev,'-termini-EPSG3413_2019-2023.shp']}; %specify consistent data sources
    end
    T_inds = []; YYYYMMDD = []; Tseas = []; Tyrs = []; Tmos = []; TX = []; TY = [];
    for l = 1:length(termfile_names)
        % load the shapefile
        % if ismember(termfiles(l).name,termfile_names)
        term = shaperead(char(termfile_names(l)));
        for k = 1:length(term)
            term_date(k) = datenum(term(k).Date,'yyyy-mm-dd');
        end
        [~,idx] = sort(term_date);
        for k = 1:length(term)
            sorted_term(k) = term(idx(k));
        end
        clear term; term = sorted_term; clear sorted_term term_date;
        %convert format of TermPicks dates
        bad_ind = [];
        for k = 1:length(term)
            %isolate the year, moth, and day from the date
            if ~isfield(term(k),'Year')
                term(k).Year = num2str(year(term(k).Date));
                term(k).Month = num2str(month(term(k).Date));
                term(k).Day = num2str(day(term(k).Date));
            else
                if k == 1 %wipe out existing data in case format is inconsistent with what I want
                    term = rmfield(term,{'Year','Month','Day'});
                end
                term(k).Year = num2str(year(term(k).Date));
                term(k).Month = num2str(month(term(k).Date));
                term(k).Day = num2str(day(term(k).Date));
            end

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

        %concatenate data from different datasets
        if l == 1
            cat_ref = 0;
        else
            cat_ref = length(T);
        end
        for k = 1:length(term)
            T(cat_ref+k).X = term(k).X; T(cat_ref+k).Y = term(k).Y;
            Tyrs(cat_ref+k,1) = str2num(term(k).Year);
            Tmos(cat_ref+k,1) = str2num(term(k).Month);
            for seas = 1:4
                if ismember(str2num(term(k).Month),seasons(seas,:))
                    Tseas(cat_ref+k,1) = seas;
                end
            end
        end
        %clear variables
        clear term YYYYMMDD datestr idx seas;
        % end
    end


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
    map_fig = figure; set(map_fig,'position',[850 50 800 600]); ax1 = axes;
    % imagesc(ax1,im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset,[],[])); axis xy equal; 
    imagesc(ax1,im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),im_subset); axis xy equal; 
    colormap(ax1,'gray'); drawnow; hold on;
    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]);
    %plot dummy lines to create a legend
    first_term = find(term_trace == 1,1,'first');
    % for l = 1:4
    %     pt(l) = plot(melmask.dated(first_term).x,melmask.dated(first_term).y,...
    %         '-','color',seas_cmap(l,:),'linewidth',1.5); hold on;
    % end

    %melange extent plotting: option 1 = plot the number of iceberg observations at each point to outline the
    %extent of the melange (instead of using melmask polygon)
    imx = single(im.x(min(xlims):max(xlims))); imy = single(im.y(min(ylims):max(ylims)));
    [Xgrid,Ygrid] = meshgrid(imx,imy); %convert the vector coordinates into matrices

    %loop through the DEMs & create a DEM mosaic
    disp('looping through the DEMs & loading elevations to a 3D matrix');
    melange_mats = dir([root_dir,site_abbrev,'/DEMs/*DEMfilled.mat']);
    for p = 1:length(melange_mats)
        DEM_name = melange_mats(p).name;
        load([root_dir,site_abbrev,'/DEMs/',DEM_name]);
        % figure; imagesc(M.DEM.x,M.DEM.y,M.DEM.z_filled); axis xy equal; hold on;
        % colormap jet;
        % plot(melmask.dated(p).x,melmask.dated(p).y,'-k','linewidth',2); hold on;
        [ZXgrid,ZYgrid] = meshgrid(M.DEM.x,M.DEM.y);

        %isolate the data in the melange mask
        DEMmask = zeros(size(ZXgrid));
        maskedz = M.DEM.z_filled.*M.mask.DEM;
        DEMmask(maskedz > zthresh) = 1;
        % figure; imagesc(M.DEM.x,M.DEM.y,DEMmask); axis xy equal; hold on;
        % colormap gray;

        %interpolate to the cropped Landsat image
        mel_data(:,:,p) = single(interp2(ZXgrid,ZYgrid,double(DEMmask),Xgrid,Ygrid)); %interpolate elevations to the Landsat 8 coordinates
        clear Z* M DEMmask maskedz in;
    end
    melange_obs = sum(mel_data,3,"omitnan"); 
    meltrans = melange_obs; meltrans(meltrans>0) = 1;
    % elev_cmap = cmocean('thermal',p); 
    % % elev_cmap = colormap(pink(p)); 
    white_end = [231,212,232]/255; purple_end = [64,0,75]/255; colorgrad = (white_end-purple_end)/(p+1);
    for k = 1:p+1
        elev_cmap(k,:) = white_end-(k-1)*colorgrad;
    end
    % elev_cmap = flipud(elev_cmap); elev_cmap(1,:) = [0 0 0]; %black
    elev_cmap(1,:) = [1 1 1]; %white

    %plot the time-averaged DEM & redraw the melange mask
    % figure1 = figure; set(gcf,'position',[50 50 1600 600]);
    ax2 = axes;
    implot = imagesc(ax2,imx,imy,melange_obs/p,'AlphaData', meltrans); axis xy equal; hold on;
    colormap(ax2,elev_cmap); set(ax2,'color','none','visible','off');
    cbar = colorbar; cbar.YLabel.String = 'Iceberg cover frequency';
    cbar.FontSize = 16; cbar.Location = "eastoutside";
    colormap(ax1,'gray');
    linkaxes([ax1 ax2]);
    % plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',2);
    clear imx imy Xgrid Ygrid DEMmask;

    %melange extent plotting: option 2 = plot faint colored boxes for the DEM binning
    % for k = min(seaward_ext):max(inland_ext)
    %     fill([MP(j).Z.transectX(k,1),MP(j).Z.transectX(k,2),MP(j).Z.transectX(k+1,2),MP(j).Z.transectX(k+1,1),MP(j).Z.transectX(k,1)],...
    %         [MP(j).Z.transectY(k,1),MP(j).Z.transectY(k,2),MP(j).Z.transectY(k+1,2),MP(j).Z.transectY(k+1,1),MP(j).Z.transectY(k,1)],...
    %         tran_cmap(k-min(seaward_ext)+1,:),'FaceAlpha',0.5,'EdgeColor','none'); hold on;
    % end

    %plot the terminus positions
    termdists = MP(j).Z.termdist; termdists(term_trace == 0) = NaN;
    [~,maxind] = max(termdists); [~,minind] = min(termdists);
    plot(melmask.dated(minind).x,melmask.dated(minind).y,'-','color',[90,174,97]/255,'linewidth',3); hold on;
    plot(melmask.dated(maxind).x,melmask.dated(maxind).y,'-','color',[90,174,97]/255,'linewidth',3); hold on;
    disp('plotted min & max centerline terminus positions in green');
    % for p = 1:size(melmask.dated,2)
    %     if term_trace(p) == 1
    %         plot(melmask.dated(p).x,melmask.dated(p).y,'-','color',...
    %             'k','linewidth',1.5); hold on;
    %         % plot(melmask.dated(p).x,melmask.dated(p).y,'-','color',...
    %         %     seas_cmap(zseas(p),:),'linewidth',1.5); hold on;
    %         drawnow;
    %     end
    % end
    % for k = 1:length(T)
    %     if Tyrs(k) >= min(zyrs) && Tyrs(k) <= max(zyrs)
    %     in = inpolygon(T(k).X,T(k).Y,melmask.uncropped.x,melmask.uncropped.y);
    %     % plot(MP(j).T.X(in),MP(j).T.Y(in),'-','color',...
    %     %     yr_cmap(Tyrs(k)-min(years)+1,:),'linewidth',1.5); hold on;
    %     plot(T(k).X(in),T(k).Y(in),'-','color',...
    %         seas_cmap(Tseas(k),:),'linewidth',1.5); hold on;
    %     clear in; drawnow;
    %     end
    % end
    %add centerline & transect locations
    plot(MP(j).V.X(min(seaward_ext):max(inland_ext)+1),MP(j).V.Y(min(seaward_ext):max(inland_ext)+1),'.-k','markersize',16,'linewidth',2); hold on;
    xlims = [min(MP(j).Z.transectX(min(seaward_ext):max(inland_ext)+1,:),[],'all'), max(MP(j).Z.transectX(min(seaward_ext):max(inland_ext)+1,:),[],'all')];
    ylims = [min(MP(j).Z.transectY(min(seaward_ext):max(inland_ext)+1,:),[],'all'), max(MP(j).Z.transectY(min(seaward_ext):max(inland_ext)+1,:),[],'all')];
    set(ax1,'xlim',xlims,'ylim',ylims);
    xticks = get(ax1,'xtick'); yticks = get(ax1,'ytick');
    set(ax1,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
    xlabel(ax1,'Easting (km)','fontsize',16); ylabel(ax1,'Northing (km)','fontsize',16);
    pos = get(gca,'position');
    if range(xlims) > 1.05*range(ylims) %short and fat map so plot the legend below
        % map_leg = legend(pt,season_names,'Location','southoutside',...
        %     'Orientation','horizontal');
        % leg_pos = get(map_leg,'position');
        % set(map_leg,'position',[leg_pos(1)+(0.5-mean([leg_pos(1) leg_pos(1)+leg_pos(3)])) 0.075 leg_pos(3) leg_pos(4)]);
        if range(xlims)/range(ylims) < 1.5
            set(ax1,'position',[pos(1) 0.25 0.9*pos(3) 0.9*pos(4)]); drawnow;
            set(ax2,'position',[pos(1) 0.25 0.9*pos(3) 0.9*pos(4)]); drawnow;
        else
            set(ax1,'position',[pos(1) 0.2 pos(3) pos(4)]); drawnow;
            set(ax2,'position',[pos(1) 0.2 pos(3) pos(4)]); drawnow;
        end
    else %tall and thin map so plot the legend on the side
        % map_leg = legend(pt,season_names,'Location','eastoutside',...
        %     'Orientation','vertical');
        % leg_pos = get(map_leg,'position');
        % if leg_pos(1)+0.05 > pos(1)-0.05+pos(3)
        %     set(map_leg,'position',[leg_pos(1)+0.05 leg_pos(2) leg_pos(3) leg_pos(4)]);
        % else
        %     set(map_leg,'position',[pos(1)-0.05+pos(3) leg_pos(2) leg_pos(3) leg_pos(4)]);
        % end
        set(ax1,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]); drawnow;
        set(ax2,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]); drawnow;
    end
    saveas(map_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-site-map.png'],'png'); %save the image
    % uiwait %advance only after figure is closed

    clear xlims ylims pt elev_cmap maskedz meltrans mel_data termdists maxind minind;
    clear berg_* im* inland_ext LCdir leg* map* melmask zmos on pos ref* seaward_ext size_classes sort_ind T T_inds term_* Tyrs Tmos *lims zyrs zdate zseas Tseas tran_cmap colorgrad term;
end
disp('Done creating site maps');


%% plot iceberg size distributions as a gif & as subplots
%distribution for a site, then adds the fall distribution (if it exists)
%with arrows showing change between seasons for small (~100 m^2) and big
%(~10^5 m^2) icebergs

%find size indices for arrows
si = find(MP(1).D.area <= 100,1,'last');
bi = find(MP(1).D.area <= 10^5,1,'last');

%create the figure
dist_fig = figure; set(dist_fig,'position',[50 850 800 600]);
loglog(MP(1).D.area,MP(geo_ind(1)).D.bergs(2,:),'-','color',seas_cmap(2,:),'linewidth',2); hold on;
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',16); grid on; 
xlabel('Surface area (m^2)','fontsize',16); ylabel('Normalized iceberg count','fontsize',16); 
drawnow;
nimages = 1;
for j = 1:length(geo_ind)
    if sum(~isnan(MP(geo_ind(j)).D.bergs(2,:))) > 0
        %plot spring
        loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(2,:),'-','color',seas_cmap(2,:),'linewidth',2); hold on;
        title(MP(geo_ind(j)).name);
        drawnow;
        frame = getframe(dist_fig);
        gif_im{nimages} = frame2im(frame); nimages = nimages+1;
        
        for p = 3:4
            %plot summer and/or spring then fall depending on data
            if sum(~isnan(MP(geo_ind(j)).D.bergs(p,:))) > 0
                loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(2,:),'-','color',seas_cmap(2,:),'linewidth',2); hold on;
                title(MP(geo_ind(j)).name);
                drawnow;
                
                %later distribution
                loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(p,:),'-','color',seas_cmap(p,:),'linewidth',2); hold on;
                drawnow;
                %             frame = getframe(dist_fig);
                %             gif_im{nimages} = frame2im(frame); nimages = nimages+1;
                
                %add arrows for small iceberg abundance change
                sp1 = [MP(1).D.area(si) MP(geo_ind(j)).D.bergs(2,si)];
                sp2 = [MP(1).D.area(si) MP(geo_ind(j)).D.bergs(p,si)];
                sdp = sp2-sp1;
                quiver(sp1(1),sp1(2),sdp(1),sdp(2),0,'linewidth',2,'color','k'); hold on;
                if sdp(2) > 0
                    plot(sp2(1),sp2(2)-0.2*sp2(2),'^k','markerfacecolor','k'); hold on;
                else
                    plot(sp2(1),sp2(2)+0.2*sp2(2),'vk','markerfacecolor','k'); hold on;
                end
                clear sp1 sp2 dsp;
                %add arrows for large iceberg abundance change
                sp1 = [MP(1).D.area(bi) MP(geo_ind(j)).D.bergs(2,bi)];
                sp2 = [MP(1).D.area(bi) MP(geo_ind(j)).D.bergs(p,bi)];
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

filename = [root_dir,'GrIS-melange_seasonal-iceberg-size-distributions_timeseries.gif']; % Specify the output file name
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
for j = 1:length(geo_ind)
        %plot spring
        subplot(sub1);
        loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(2,:),'-','color',seas_cmap(2,:),'linewidth',2); hold on;
        drawnow;
        
        %plot summer
        subplot(sub2);
        loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(3,:),'-','color',seas_cmap(3,:),'linewidth',2); hold on;
        
        %plot fall
        subplot(sub3);
        loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(4,:),'-','color',seas_cmap(4,:),'linewidth',2); hold on;
        drawnow;

end
subplot(sub1);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',16); grid on; 
text(10^5.8,10^-1,'spring','fontsize',16)
xlabel('Surface area (m^2)','fontsize',16); ylabel('Normalized iceberg count','fontsize',16);
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
subplot(sub2);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'yticklabel',[],'fontsize',16); grid on;
text(10^5.8,10^-1,'summer','fontsize',16)
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
subplot(sub3);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'yticklabel',[],'fontsize',16); grid on; 
text(10^5.8,10^-1,'fall','fontsize',16)
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
drawnow;
saveas(dist_fig,[root_dir,'GrIS-melange_iceberg-distribution_loglog.png'],'png'); %save the plot
close all;

%% 

%Size distribution parameters
nthresh = 1e-6; % set small number bin cutoff (n1 must be greater than this value)
zthresh = 3; %set small size bin cutoff (freeboard must exceed this value)
rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)
ARcomp.best.autoALL = 2; % iceberg aspect ratio 
vthresh = (1/4)*pi*((rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL.*zthresh).^2; %don't include this bin size or smaller in curve fitting
dplawthresh = 10^5; % upper bound on the intercept for the dummy powerlaw
norm_type = 2; % toggle between L2, max, and log norm using 2, Inf, and 'log'
normalize_exp = 1.5; % Increase to weight residuals towards end of the curve, minimum = 1

%create a similar seasonal profile figure but with a top row containing
%profiles near the terminus and a bottom row containing seaward profiles
subdist_fig = figure; set(subdist_fig,'position',[50 50 1200 1200]);
%subplots for inland sample location
subia = subplot(2,3,1); subib = subplot(2,3,2); subic = subplot(2,3,3); 
%subplots for seaward sample location
subsa = subplot(2,3,4); subsb = subplot(2,3,5); subsc = subplot(2,3,6); 
%letters for the plots
alphabet = ['a','b','c','d','e','f'];
for j = 1:length(geo_ind)
    seasons = MP(geo_ind(j)).D.months; site_abbrev = MP(geo_ind(j)).name;
    % disp(sitenames(geo_ind(j),:)); 
    
    %isolate the dates for data sorting
    for p = 1:length(MP(geo_ind(j)).Z.date)
        term_trace(p) = MP(geo_ind(j)).Z.termflag(p);
        % zdate(p) = convert_to_decimaldate(char(MP(geo_ind(j)).Z.date(p)));
        % datest = datetime(MP(geo_ind(j)).Z.date{p},'InputFormat','yyyyMMdd');
        % zyrs(p) = year(datest); zmos(p) = month(datest);
        % clear datest;
    end

    %loop through subset size distributions & find the first and last
    %non-NaN columns to identify melange extent
    Dsubs = dir([root_dir,MP(geo_ind(j)).name,'/',MP(geo_ind(j)).name,'*-iceberg-distribution-subsets.csv']);
    size_classes = [];
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(geo_ind(j)).name,'/',Dsubs(p).name],"VariableNamingRule","preserve");
        name_split = split(Dsubs(p).name,'-',2);
        berg_dates(p) = datetime(name_split{2},'Inputformat','yyyyMMdd'); clear name_split;
        mos(p) = month(berg_dates(p));

        %identify observational limits along the centerline
        berg_nos = table2array(D(:,3:end)); berg_nos(berg_nos==0) = NaN;
        size_classes(p,:) = sum(~isnan(berg_nos),1);
        seaward_ext(p) = find(size_classes(p,:)>0,1,'first');
        inland_ext(p) = find(size_classes(p,:)>0,1,'last')+1; inland_idx(p) = inland_ext(p);
        clear berg_nos;

        %fill-in inland boundary reference based on terminus observations
        %from images if DEM did not overlap the terminus
        if term_trace(p) == 0
            if ~isnan(MP(geo_ind(j)).Z.termdist(p))
                dists = sqrt((MP(geo_ind(j)).Z.termX(p)-MP(geo_ind(j)).V.X).^2 + (MP(geo_ind(j)).Z.termY(p)-MP(geo_ind(j)).V.Y).^2);
                [~,ia] = sort(dists);
                inland_idx(p) = max(ia(1:2)); term_trace(p) = 1; %use terminus traces from imagery
                clear dists ia;
            else
                inland_idx(p) = NaN;
                disp(['Need terminus data for ',MP(geo_ind(j)).name,' ',char(MP(geo_ind(j)).Z.date(p))])
            end
        end

    end

    %now analyze the ize distributions
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(geo_ind(j)).name,'/',Dsubs(p).name],"VariableNamingRule","preserve");

        %extract size distribution profiles
        if p == 1 %initialize the data compilation
            v1 = table2array(D(:,1)); dv1 = table2array(D(:,2));
        end
        %first 2 columns are area & bin width, so 3+ are data from points along the centerline
        n1 = NaN(length(v1),max(inland_idx)-1);
        if ~isnan(inland_idx(p))
            berg_nos = fliplr(table2array(D(:,3:(inland_idx(p)-1)+2)));
            berg_nos(berg_nos==0) = NaN;
            n1(:,1:size(berg_nos,2)) = berg_nos./dv1;
        end
        clear D;
        %stack distributions in a structure
        binned_n1_stack(:,:,p) = n1;
        clear n1 berg_nos;
    end

    %assign seaward & inland sampling limits based on the selected sampling
    %strategy (sampling = ['fixed','dated']
    inland_meanidx = round(mean(inland_ext(term_trace==1)));
    seaward_meanidx = round(mean(seaward_ext));
    seaward_idx = inland_meanidx-seaward_meanidx;

    %create seasonal average distributions for each fjord bin
    for k = 1:4
        % disp([char(season_names{k}),': parameters starting from the terminus & moving seaward...']);

        %stack the profiles
        n1_seasmean(:,:,k) = nanmean(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),3);
        % n1_seasmax(:,:,k) = (nanmean(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),3)+std(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),0,3,'omitnan'));
        % n1_seasmin(:,:,k) = (nanmean(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),3)-std(binned_n1_stack(:,:,ismember(mos,seasons(k,:))==1),0,3,'omitnan'));

        for l = [1,seaward_idx]
            n1 = n1_seasmean(:,l,k);
            if l == 1
                n1_inlandseas(:,j,k) = n1;
            else
                n1_seawardseas(:,j,k) = n1;
            end
            clear n1;
        end

        % %plot the seasonal profiles at each location
        % if k == 2; subplot(subia); elseif k == 3; subplot(subib); else; subplot(subic); end
        % loglog(v1,n1_inlandseas(:,j,k),'-','LineWidth',1,'Color',seas_cmap(k,:)); hold on;
        % if k == 2; subplot(subsa); elseif k == 3; subplot(subsb); else; subplot(subsc); end
        % loglog(v1,n1_seawardseas(:,j,k),'-','LineWidth',1,'Color',seas_cmap(k,:)); hold on;

    end
    %concatenate the profiles for each season so that they can be plotted
    %as filled polygons instead of individual lines
    %seaward
    n1wi(:,j,1) = n1_seawardseas(:,j,1); n1sp(:,j,1) = n1_seawardseas(:,j,2); 
    n1su(:,j,1) = n1_seawardseas(:,j,3); n1fa(:,j,1) = n1_seawardseas(:,j,4);
    %inland
    n1wi(:,j,2) = n1_inlandseas(:,j,1); n1sp(:,j,2) = n1_inlandseas(:,j,2); 
    n1su(:,j,2) = n1_inlandseas(:,j,3); n1fa(:,j,2) = n1_inlandseas(:,j,4);

    clear berg_dates mos berg_nos size_classes seaward_* inland_* binned_n1_stack n1_seasmean term_trace;
end
%plot a filled polygon to bound the observed size distributions
log_v = log(v1); 
subplot(subia); %inland, spring
log_n = log(n1sp); maxn = max(log_n(:,:,2),[],1); minn = min(log_n(:,:,2),[],1); 
fill([log_v, flipud(log_v)],[maxn; flipud(minn)],...
    seas_cmap(2,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
clear log_n maxn minn;
subplot(subsa); %seaward, spring
maxn = max(log_n(:,:,1),[],1); minn = min(log_n(:,:,1),[],1); 
fill([log_v, flipud(log_v)],[maxn; flipud(minn)],...
    seas_cmap(2,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
clear log_n maxn minn;
subplot(subib); %inland, summer
log_n = log(n1su); maxn = max(log_n(:,:,2),[],1); minn = min(log_n(:,:,2),[],1); 
fill([log_v, flipud(log_v)],[maxn; flipud(minn)],...
    seas_cmap(3,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
clear log_n maxn minn;
subplot(subsb); %seaward, summer
maxn = max(log_n(:,:,1),[],1); minn = min(log_n(:,:,1),[],1); 
fill([log_v, flipud(log_v)],[maxn; flipud(minn)],...
    seas_cmap(3,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
clear log_n maxn minn;
subplot(subic); %inland, fall
log_n = log(n1fa); maxn = max(log_n(:,:,2),[],1); minn = min(log_n(:,:,2),[],1); 
fill([log_v, flipud(log_v)],[maxn; flipud(minn)],...
    seas_cmap(4,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
clear log_n maxn minn;
subplot(subsc); %seaward, fall
maxn = max(log_n(:,:,1),[],1); minn = min(log_n(:,:,1),[],1); 
fill([log_v, flipud(log_v)],[maxn; flipud(minn)],...
    seas_cmap(4,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
clear log_n maxn minn;
drawnow;
%apply the EBC fit function to the average of all seasonal profiles in each
%location & plot the result as a dashed line
for l = 1:2
    if l == 1
        disp('plotting profiles for near-terminus sample location')
    else
        disp('plotting profiles for seaward sample location')
    end

    for k = 2:4
        if l == 1
            n1 = nanmean(n1_inlandseas(:,:,k),2);
        else
            n1 = nanmean(n1_seawardseas(:,:,k),2);
        end

        %filter to decide if there are enough points to fit lines
        v = v1(n1>nthresh); n = n1(n1>nthresh); dv = dv1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
        n = n(v>vthresh); dv = dv(v>vthresh); v = v(v>vthresh); %remove smallest size classes
        v = v(~isnan(n)); dv = dv(~isnan(n)); n = n(~isnan(n)); % Remove NaNs

        % Fit E-BC fragmentation model to the data and plot
        if sum(~isnan(v)) > 12
            %enter unfiltered data into model (filtering happens there)
            tempfig = figure;
            [alpha,c1,c2,c3,c4,~,~] = EBC_fragmentation_curve(v1, n1, dv1, norm_type, vthresh, nthresh, 1, normalize_exp); % fit Eq. 1
            n1_mod = EBC_model([c1,c2,alpha,c3,c4],v1); % grab model
            close(tempfig);

            %plot the size distribution fit
            figure(subdist_fig);
            if l == 1
                if k == 2; subplot(subia); elseif k == 3; subplot(subib); else; subplot(subic); end
                
                %plot the data
                loglog([v1(zthresh+1),v1(zthresh+1)],[10^-12 10^3],'-','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                loglog(v1,n1,'.-','LineWidth',2,'Color','k'); hold on;
                loglog(v1,n1_mod,'--','LineWidth',2,'Color','k'); hold on;

                %format the subplot
                set(gca,'ylim',[10^-10 10^3],'xlim',[10^1 10^7],'xtick',[10^2, 10^4, 10^6],'fontsize',16); grid on;
                text(10^5,10^2,[alphabet(k*l-1),') ',char(season_names(k))],'fontsize',16)
                if k == 2
                    ylabel('Normalized iceberg count','fontsize',16);
                end

            else
                if k == 2; subplot(subsa); elseif k == 3; subplot(subsb); else; subplot(subsc); end
                
                %plot the data
                loglog([v1(zthresh+1),v1(zthresh+1)],[10^-12 10^3],'-','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                loglog(v1,n1,'.-','LineWidth',2,'Color','k'); hold on;
                loglog(v1,n1_mod,'--','LineWidth',2,'Color','k'); hold on;

                %format the subplot
                set(gca,'ylim',[10^-10 10^3],'xlim',[10^1 10^7],'xtick',[10^2, 10^4, 10^6],'fontsize',16); grid on;
                text(10^5,10^2,[alphabet((k-1)+3),') ',char(season_names(k))],'fontsize',16)
                xlabel('Surface area (m^2)','fontsize',16);
                if k == 2
                    ylabel('Normalized iceberg count','fontsize',16);
                end
            end

            clear n1_mod v dv n n1 alpha c1 c2 c3 c4;
        end
        drawnow;
    end
    disp('moving on')
end
saveas(subdist_fig,[root_dir,'GrIS-melange_iceberg-distribution-subsets_loglog.png'],'png'); %save the plot

%% load TermPicks timeseries for each glacier
% cd(root_dir);
% load([root_dir,'GrIS-melange-characteristics.mat']);
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
%     save([root_dir,'GrIS-melange-characteristics.mat'],'MP','-v7.3');
% 
%     %clear variables
%     clear term YYYYMMDD datestr idx site_abbrev;
% end



%% fit iceberg size distribution curves to average of normalized, binned, size distributions





