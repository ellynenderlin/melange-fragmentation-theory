%%% Compile Greenland melange characteristics into a single matfile

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

%Thickness parameters:
zcutoff = zthresh; %elevation threshold below which to ignore iceberge (m)
rho_i = 900; %ice density (kg/m^3)
rho_w = 1026; %water density (kg/m^3)
Hcutoff = round((rho_w/(rho_w-rho_i))*zcutoff); %H threshold for figure naming

%Size distribution parameters
nthresh = 1e-6; % set small number bin cutoff (n1 must be greater than this value)
% zthresh = 3; %set small size bin cutoff (freeboard must exceed this value)
rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)
ARcomp.best.autoALL = 2; % iceberg aspect ratio 
vthresh = (1/4)*pi*((rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL.*zthresh).^2; %don't include this bin size or smaller in curve fitting
dplawthresh = 10^5; % upper bound on the intercept for the dummy powerlaw
norm_type = 2; % toggle between L2, max, and log norm using 2, Inf, and 'log'
normalize_exp = 1.5; % Increase to weight residuals towards end of the curve, minimum = 1

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
geo_names = [{'Ullip'},{'Nuussuup'},{'Nunatakassaap'},...
    {'Illullip'},{'Upernavik North'},{'Salliarutsip'},{'Umiammakku'},...
    {'Kangilliup'},{'Sermeq Kujalleq'},{'Helheim'},{'Midgard'},...
    {'Magga Dan'},{'Daugaard Jensen'},{'Zachariae Isstrom'}];
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
                disp(['...restarting on site #',num2str(site_start),' (',sitenames(site_start,:),')']);
            else
                disp('data reloaded and dataset is fully processed');
            end
        case '2) No: start fresh'
            site_start = 1;
            MP = struct;
    end
else
    site_start = 1;
    MP = struct;
end
disp('Move on to the next subsection in order to compile data & generate figures.')

%% loop through the folders & extract info

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

    %Extract width-averaged elevations from the transects
    %find the NaNs in the coordinate pairs to identify each transect
    coords = table2array(T(:,1:2));
    nan_inds = find(isnan(coords(:,1))==1);

    %extract the data
    tran_ind = 1; k = 1; zprofs = [];
    while k < size(T,1)
        if k == 1
            %not needed for buttressing estimates but can be used for plotting melange bins: 
            %extract full coordinates for transect
            MP(j).Z.transectX(tran_ind,1:2) = coords([1,nan_inds(k)-1],2);
            MP(j).Z.transectY(tran_ind,1:2) = coords([1,nan_inds(k)-1],1);

            % %not needed for buttressing estimates: extract the centroid coordinates
            % MP(j).Z.centerX(tran_ind,1) = nanmean(table2array(T(1:nan_inds(k)-1,2)));
            % MP(j).Z.centerY(tran_ind,1) = nanmean(table2array(T(1:nan_inds(k)-1,1)));

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


    %Not needed for buttressing estimates but can use to grab velocities corresponding to each DEM date:
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
            short_dts = find(V.days_dt>vdtmin & V.days_dt<vdtmax); %get rid of all velocities with coarse temporal resolution
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

    %standardize the x-limits on the plots
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


%% extract melange attributes, estimate buttressing, & make overview plots
close all; drawnow;

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

% %If adding new dates, profiles need to be removed from the structure
% %because the inland reference may have changed
% %(will throw errors when concatenating)
% for j = 1:length(MP)
%     MP(j).Z = rmfield(MP(j).Z,{'dist','Hseas'});
%     MP(j).V = rmfield(MP(j).V,{'dist','Vseas','dVdx'});
% end

%create dummy vectors to hold all buttressing data for histogram plotting
BM_annual = []; BM_character = []; BA_annual = []; BA_character = [];
site_naming = []; PS_X = []; PS_Y = []; sp_stats = [];

disp(['Creating profile plots ignoring icebergs thinner than ',num2str(Hcutoff),'m & using ',vfilter,' speeds w/ dt ~',num2str(round([vdtmin,vdtmax])),' days']);
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
    Havg(term_trace==0,:) = NaN; H_seas = NaN(max(inland_idx)-1,4,length(years)); pack_seas = NaN(max(inland_idx)-1,4,length(years)); 
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
        byrs(p) = year(berg_dates(p)); bmos(p) = month(berg_dates(p));

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
    %dummy matrices to hold seasonal average distributions for each year
    bergdist_inland_seas = NaN(4,size(bergdist_inland,2),length(years));
    bergdist_seaward_seas = NaN(4,size(bergdist_inland,2),length(years));
    for p = 1:length(years)
        yr_idx = find(byrs == years(p));
        if ~isempty(yr_idx)

            %isolate the distributions for that year
            for k = 1:length(yr_idx)
                bi_temp(k,:) = bergdist_inland_norm(yr_idx(k),:);
                bs_temp(k,:) = bergdist_seaward_setdx_norm(yr_idx(k),:); %could change to bergdist_seaward_end_norm for moving end
            end

            %calculate seasonal average
            mos_yr = bmos(yr_idx);
            for k = 1:4
                if ~isempty(ismember(mos_yr,seasons(k,:)))
                    bergdist_inland_seas(k,:,p) = nanmean(bi_temp(ismember(mos_yr,seasons(k,:))==1,:),1); %inland
                    bergdist_sewaward_seas(k,:,p) = nanmean(bs_temp(ismember(mos_yr,seasons(k,:))==1,:),1); %seaward
                end
            end

            clear mos_yr b*_temp;
        end
        clear yr_idx
    end
    %average the seasonal average datasets
    bergdist_seas = NaN(4,size(bergdist_inland,2),2);
    for k = 1:4
        if ~isempty(ismember(bmos,seasons(k,:)))
            bergdist_seas(k,:,1) = nanmean(bergdist_inland_seas(k,:,:),3); %inland
            bergdist_seas(k,:,2) = nanmean(bergdist_sewaward_seas(k,:,:),3); %seaward
        end
    end
    disp('extracted distributions');

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
            binned_n1_stack(:,:,p) = n1; %rows = counts, columns = bin wrt terminus, p = date
            clear n1 berg_nos;
        end
        %create seasonal average distributions for each fjord bin
        n1wi = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
        n1sp = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
        n1su = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
        n1fa = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
        for p = 1:length(years)
            yr_idx = find(byrs == years(p));
            if ~isempty(yr_idx)

                %isolate the distributions for that year
                for k = 1:length(yr_idx)
                    n1_temp(:,:,k) = binned_n1_stack(:,:,yr_idx(k));
                end

                %calculate seasonal averages across the year
                mos_yr = bmos(yr_idx);
                if ~isempty(ismember(mos_yr,seasons(1,:))) %winter
                    n1wi(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(1,:))==1),3); 
                end
                if ~isempty(ismember(mos_yr,seasons(2,:))) %spring
                    n1sp(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(2,:))==1),3); 
                end
                if ~isempty(ismember(mos_yr,seasons(3,:))) %summer
                    n1su(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(3,:))==1),3); 
                end
                if ~isempty(ismember(mos_yr,seasons(4,:))) %fall
                    n1fa(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(4,:))==1),3); 
                end

                clear mos_yr n1_temp;
            end
            clear yr_idx
        end
        %calculate seasonal averages across all years & fit curves
        n1_seasmean(:,:,1) = nanmean(n1wi,3);
        n1_seasmax(:,:,1) = (nanmean(n1wi,3)+std(n1wi,0,3,'omitnan')); n1_seasmin(:,:,1) = (nanmean(n1wi,3)-std(n1wi,0,3,'omitnan'));
        n1_seasmean(:,:,2) = nanmean(n1sp,3);
        n1_seasmax(:,:,2) = (nanmean(n1sp,3)+std(n1sp,0,3,'omitnan')); n1_seasmin(:,:,2) = (nanmean(n1sp,3)-std(n1sp,0,3,'omitnan'));
        n1_seasmean(:,:,3) = nanmean(n1su,3);
        n1_seasmax(:,:,3) = (nanmean(n1su,3)+std(n1su,0,3,'omitnan')); n1_seasmin(:,:,3) = (nanmean(n1su,3)-std(n1su,0,3,'omitnan'));
        n1_seasmean(:,:,4) = nanmean(n1fa,3);
        n1_seasmax(:,:,4) = (nanmean(n1fa,3)+std(n1fa,0,3,'omitnan')); n1_seasmin(:,:,4) = (nanmean(n1fa,3)-std(n1fa,0,3,'omitnan'));
        clear n1wi n1sp n1su n1fa;
        for k = 1:4
            disp([char(season_names{k}),': parameters starting from the terminus & moving seaward...']);
            
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

  
    %extract thickness profiles relative to moving terminus
    for k = 1:4
        Hmean = nanmean(H_seas(:,k,:),3);
        Hmax = (nanmean(H_seas(:,k,:),3)+std(H_seas(:,k,:),0,3,'omitnan')); 
        Hmin = (nanmean(H_seas(:,k,:),3)-std(H_seas(:,k,:),0,3,'omitnan')); 
        if sum(~isnan(Hmean)) > 0
            MP(j).Z.Hseas(k,:) = Hmean';
        else
           MP(j).Z.Hseas(k,:) = NaN(1,size(H_seas,1));
        end
        clear Hmean Hmax* Hmin*;
    end
        
    %extract speed profiles relative to moving terminus
    for k = 1:4
        v_mean = nanmean(vel_seas(:,k,:),3)./365;
        vmax = (nanmean(vel_seas(:,k,:),3)+std(vel_seas(:,k,:),0,3,'omitnan'))./365; 
        vmin = (nanmean(vel_seas(:,k,:),3)-std(vel_seas(:,k,:),0,3,'omitnan'))./365; 
        
        %mask out values near zero because they are so much slower than the
        %melange velocities that they must be bad sea ice tracking
        bottomout = find(v_mean < 0.5,1,'first');
        v_mean(bottomout:end) = NaN; vmax(bottomout:end) = NaN; vmin(bottomout:end) = NaN; 
        if sum(~isnan(v_mean)) ~= 0
            MP(j).V.Vseas(k,:) = v_mean';
        else
            MP(j).V.Vseas(k,:) = NaN(1,size(vel_seas,1));
        end
        clear v_mean vmax* vmin*;
    end
        
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
                %if no velocity for the first melange point, look one
                %further down-fjord
                if isnan(MP(j).B.dVdx(1,k,p)); MP(j).B.dVdx(1,k,p) = (diff([vel_seas(1+ref_adjust,k,p),vel_seas(3+ref_adjust,k,p)])./diff([vdist(1+ref_adjust),vdist(3+ref_adjust)]))/365; end
                press = 0.5*rho_i*(1-(rho_i/rho_w))*9.81*MP(j).B.Ho(zcutoff+1,k,p);
                MP(j).B.butt_Meng(zcutoff+1,k,p) = press*MP(j).B.packing(zcutoff+1,k,p)*MP(j).B.Ho(zcutoff+1,k,p);
                MP(j).B.butt_Amundson(zcutoff+1,k,p) = (-2*(MP(j).B.Ho(zcutoff+1,k,p)*press*MP(j).B.dVdx(1,k,p))/((MP(j).B.dVdx(1,k,p)/0.3)+MP(j).B.dVdx(1,k,p)))+press*MP(j).B.Ho(zcutoff+1,k,p);
                clear press dvdist_temp;
            end
        end
    end
    MP(j).B.butt_Meng(MP(j).B.butt_Meng==0) = NaN; MP(j).B.butt_Amundson(MP(j).B.butt_Amundson==0) = NaN; 

    %display characteristics
    disp([site_abbrev,' Melange Characteristics:'])
    disp(' Near-Terminus Thickness: m')
    disp(['    winter = ',num2str(round(nanmean(MP(j).B.Ho(zcutoff+1,1,:)),1)),', spring = ',num2str(round(nanmean(MP(j).B.Ho(zcutoff+1,2,:)),1)),', summer = ',num2str(round(nanmean(MP(j).B.Ho(zcutoff+1,3,:)),1)),', fall = ',num2str(round(nanmean(MP(j).B.Ho(zcutoff+1,4,:)),1))]);
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
    disp(['    winter = ',num2str(round(nanmean(MP(j).B.packing(zcutoff+1,1,:)),2)),', spring = ',num2str(round(nanmean(MP(j).B.packing(zcutoff+1,2,:)),2)),', summer = ',num2str(round(nanmean(MP(j).B.packing(zcutoff+1,3,:)),2)),', fall = ',num2str(round(nanmean(MP(j).B.packing(zcutoff+1,4,:)),2))]);
    disp(' Near-Terminus Bergy Bit Overestimation: % Area')
    disp(['    winter = ',num2str(round(res_A(1,1),2)),', spring = ',num2str(round(res_A(2,1),2)),', summer = ',num2str(round(res_A(3,1),2)),', fall = ',num2str(round(res_A(4,1),2))]);
    disp(' Extension-Compression Strain Rate Change Location: km')
    disp(['    winter = ',num2str(round(sr_dsign_pt(1)/1000,2)),', spring = ',num2str(round(sr_dsign_pt(2)/1000,2)),', summer = ',num2str(round(sr_dsign_pt(3)/1000,2)),', fall = ',num2str(round(sr_dsign_pt(4)/1000,2))]);
    disp(' Change in Bergy Bit Misprediction Along Flow: % Area')
    disp(['    winter = ',num2str(round(res_A(1,min(res_Alim))-res_A(1,1),2)),', spring = ',num2str(round(res_A(2,min(res_Alim))-res_A(2,1),2)),', summer = ',num2str(round(res_A(3,min(res_Alim))-res_A(3,1),2)),', fall = ',num2str(round(res_A(4,min(res_Alim))-res_A(4,1),2))]);
    disp(' ');
    disp(' Buttressing estimated using near-terminus observations:');
    disp('   Meng Eqn (thickness- & packing-based): x10^6 N/m');
    disp(['    winter = ',num2str(round(nanmean(MP(j).B.butt_Meng(zcutoff+1,1,:))/10^6,2)),', spring = ',num2str(round(nanmean(MP(j).B.butt_Meng(zcutoff+1,2,:))/10^6,2)),', summer = ',num2str(round(nanmean(MP(j).B.butt_Meng(zcutoff+1,3,:))/10^6,2)),', fall = ',num2str(round(nanmean(MP(j).B.butt_Meng(zcutoff+1,4,:))/10^6,2))]);
    disp('   Amundson Eqn (thickness- & strainrate-based): x10^6 N/m');
    disp(['    winter = ',num2str(round(nanmean(MP(j).B.butt_Amundson(zcutoff+1,1,:))/10^6,2)),', spring = ',num2str(round(nanmean(MP(j).B.butt_Amundson(zcutoff+1,2,:))/10^6,2)),', summer = ',num2str(round(nanmean(MP(j).B.butt_Amundson(zcutoff+1,3,:))/10^6,2)),', fall = ',num2str(round(nanmean(MP(j).B.butt_Amundson(zcutoff+1,4,:))/10^6,2))]);
    disp(' ');
    BM_annual = [BM_annual; squeeze(MP(j).B.butt_Meng(zcutoff+1,1,:))./10^6, squeeze(MP(j).B.butt_Meng(zcutoff+1,2,:))./10^6, squeeze(MP(j).B.butt_Meng(zcutoff+1,3,:))./10^6, squeeze(MP(j).B.butt_Meng(zcutoff+1,4,:))./10^6];
    BM_character = [BM_character; nanmean(MP(j).B.butt_Meng(zcutoff+1,1,:))/10^6, nanmean(MP(j).B.butt_Meng(zcutoff+1,2,:))/10^6, nanmean(MP(j).B.butt_Meng(zcutoff+1,3,:))/10^6, nanmean(MP(j).B.butt_Meng(zcutoff+1,4,:))/10^6];
    BA_annual = [BA_annual; squeeze(MP(j).B.butt_Amundson(zcutoff+1,1,:))./10^6, squeeze(MP(j).B.butt_Amundson(zcutoff+1,2,:))./10^6, squeeze(MP(j).B.butt_Amundson(zcutoff+1,3,:))./10^6, squeeze(MP(j).B.butt_Amundson(zcutoff+1,4,:))./10^6];
    BA_character = [BA_character; nanmean(MP(j).B.butt_Amundson(zcutoff+1,1,:))/10^6, nanmean(MP(j).B.butt_Amundson(zcutoff+1,2,:))/10^6, nanmean(MP(j).B.butt_Amundson(zcutoff+1,3,:))/10^6, nanmean(MP(j).B.butt_Amundson(zcutoff+1,4,:))/10^6];

    %combine site names, coordinates, and spring thickness+buttressing into a matrix
    %for exporting as a site summary CSV
    if ~contains(MP(j).name,'KBG')
        site_naming = [site_naming; geo_order(find(geo_ind==j)), geo_names(find(geo_ind==j))];
        PS_X = [PS_X; nanmean(MP(j).Z.termX)]; PS_Y = [PS_Y; nanmean(MP(j).Z.termY)];
        if sum(~isnan(MP(j).B.butt_Amundson(zcutoff+1,2,:))) > 0
            sp_stats = [sp_stats; nanmean(MP(j).B.Ho(zcutoff+1,2,:)),...
                nanmean(cat(3,MP(j).B.butt_Meng(zcutoff+1,2,~isnan(MP(j).B.butt_Amundson(zcutoff+1,2,:))),MP(j).B.butt_Amundson(zcutoff+1,2,~isnan(MP(j).B.butt_Amundson(zcutoff+1,2,:)))),"all")/10^6];
        else
            sp_stats = [sp_stats; nanmean(MP(j).B.Ho(zcutoff+1,2,:)),...
                nanmean(MP(j).B.butt_Meng(zcutoff+1,2,:))/10^6];
        end
    end

    %export the data to tables
    yrs_temp = []; seas_temp = [];
    %thickness and packing density
    H_temp = []; pack_temp = []; V_temp = []; B_temp = [];
    Hdist(1) = 0; %changed from 0m to 1m earlier for log plotting
    for p = 1:length(years)
        for k = 1:4
            if sum(~isnan(H_seas(:,k,p))) > 0
                yrs_temp = [yrs_temp; years(p)];
                seas_temp = [seas_temp; season_names(k)];
                H_temp = [H_temp; H_seas(:,k,p)'];
                pack_temp = [pack_temp; 100*pack_seas(:,k,p)'];
            end
        end
    end
    TH = [array2table(yrs_temp), array2table(seas_temp), array2table(H_temp)];
    for l = 1:size(H_temp,2)
        bin_name(l) = {['Thickness (m): ',num2str(Hdist(l)),'m']};
    end
    column_names = ["Year","Season", bin_name]; TH.Properties.VariableNames = column_names;
    writetable(TH,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-thickness-profiles.csv']);
    TP = [array2table(yrs_temp), array2table(seas_temp), array2table(pack_temp)];
    for l = 1:size(pack_temp,2)
        bin_name(l) = {['Packing density (%): ',num2str(Hdist(l)),'m']};
    end
    column_names = ["Year","Season", bin_name]; 
    TP.Properties.VariableNames = column_names;
    writetable(TP,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-packing-profiles.csv']);
    clear H_temp pack_temp bin_name;
    %speeds
    for p = 1:length(years)
        for k = 1:4
            if sum(~isnan(H_seas(:,k,p))) > 0 %will write a NAN profile for missing speed data for a DEM observation season
                V_temp = [V_temp; vel_seas(:,k,p)'./365];
            end
        end
    end
    TV = [array2table(yrs_temp), array2table(seas_temp), array2table(V_temp)];
    for l = 1:size(V_temp,2)
        bin_name(l) = {['Speed (m/d): ',num2str(vdist(l)),'m']};
    end
    column_names = ["Year","Season", bin_name]; TV.Properties.VariableNames = column_names;
    writetable(TV,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-speed-profiles.csv']);
    clear bin_name V_temp;
    %buttressing variables and estimates
    for p = 1:length(years)
        for k = 1:4
            if sum(~isnan(H_seas(:,k,p))) > 0 %will write a NAN profile for missing speed data for a DEM observation season
                B_temp = [B_temp; round(MP(j).B.Ho(zcutoff+1,k,p),2),round(MP(j).B.packing(zcutoff+1,k,p),4),round(MP(j).B.butt_Meng(zcutoff+1,k,p),0),...
                    round(MP(j).B.dVdx(1,k,p),7),round(MP(j).B.butt_Amundson(zcutoff+1,k,p),0)];
            end
        end
    end
    TB = [array2table(yrs_temp), array2table(seas_temp), array2table(B_temp)];
    column_names = ["Year","Season","Thickness (m)","Packing (fraction)","Buttressing-Meng (N/m)",...
        "Strainrate (1/d)","Buttressing-Amundson (N/m)"]; TB.Properties.VariableNames = column_names;
    writetable(TB,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-buttressing.csv']);
    clear B_temp yrs_temp seas_temp;
    clear TH TP TV TB;
    disp('Profiles and buttressing time series saved as CSVs')

    %save the matfile
    save([root_dir,'GrIS-melange-characteristics.mat'],'MP','-v7.3');
    
    %refresh
    clear berg_* bergdist* berg_normdist* C centerline* D Dsubs *idx seaward_* inland_* term_* tran_* size_classes Zfilt H_* Havg vel_* vels* v_mean w zdate berg_mo bins bin_no z_* pos pz pv *dist *yrs *mos seas_leg packing;
    clear leg_* ref_* Tdate* xlims Tslopes res* Tres size_plslope bin* column_names dv1 v1 n1* *_seas ax ax2 max_xlim width_prof Vterm dVdxterm subv subz*;
    clear ax1 ax2;

end

%export overview data to a CSV (used for plotting sitemap in QGIS)
TS = [array2table(site_naming(:,2)), array2table(PS_X), array2table(PS_Y), array2table(sp_stats)];
column_names = ["Site Name","X (m)", "Y (m)","Melange thickness (m)","Melange buttressing (x10^6 N/m)"]; TS.Properties.VariableNames = column_names;
writetable(TS,[root_dir,'GrIS-melange-sites.csv']);
clear TS;


