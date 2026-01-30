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

%load the melange characteristic dataset created by compile_melange_characteristics.m
load([root_dir,'GrIS-melange-characteristics.mat']);
            


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
        for k = 1:4
            if ismember(zmos(p),seasons(k,:))==1
                zseas(p) = k;
            end
        end
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
        for k = 1:4
            if ismember(Tmos(p),seasons(k,:))==1
                Tseas(p) = k;
            end
        end
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
            'color',seas_cmap(Tseas(p),:),'linewidth',2,'markersize',10); hold on;
    end
    for p = 1:length(MP(j).Z.date)
        if term_trace(p) == 0 && isnan(MP(j).Z.termdist(p))
            disp(['Still need terminus data for ',MP(j).name,' ',char(MP(j).Z.date(p))]);
        else
            plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).Z.termdist(p),zdate(p),'s',...
                'color',seas_cmap(zseas(p),:),'linewidth',1,'markerfacecolor','none','markersize',10); hold on;
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

    %create characteristic seasonal curves for terminus position
    MP(j).T.termanom = NaN(size(MP(j).T.termdist));
    for p = 1:length(years)
        yr_idx = find(Tyrs == years(p));
        if ~isempty(yr_idx)
            %identify the months in the specificyear
            mos_yr = Tmos(yr_idx);

            %calculate seasonal medians
            for k = 1:4
                Tdist_seas(k,p) = nanmedian(MP(j).T.termdist(yr_idx(ismember(mos_yr,seasons(k,:))==1)));
            end

            %calculate the mean of the seasonal positions
            Tmean = nanmean(Tdist_seas(:,p));

            %subtract seasonal medians from the annual median so that
            %positive anomalies indicate terminus advance (farther from
            %most retreated terminus) & negative anomalies indicate retreat
            MP(j).T.termanom(yr_idx) = Tmean - MP(j).T.termdist(yr_idx);
            Tdist_seasanom(:,p) = Tmean - Tdist_seas(:,p);
            
            clear Tmean mos_y;
        end
        clear yr_idx;
    end

    %calculate monthly anomalies
    MP(j).T.termanom_mo = NaN(12,1); MP(j).T.termanom_morange = NaN(12,2); 
    for k = 1:12
        MP(j).T.termanom_mo(k,1) = nanmean(MP(j).T.termanom(ismember(Tmos,k)==1));
        MP(j).T.termanom_morange(k,:) = prctile(MP(j).T.termanom(ismember(Tmos,k)==1),[5,95],2);
    end

    clear C centerline_dist dts inland_idx seaward_idx mean_prof Tdate* term_* Tmos tran_dist Tyrs zdate* zmos zyrs;
    clear Tdist* zseas Tseas;
end
%resave the data (as needed)
save([root_dir,'GrIS-melange-characteristics.mat'],'MP','-v7.3');

%plot all the characteristic monthly seasonal terminus anomaly timeseries
seasterm_fig = figure; set(seasterm_fig,'position',[50 50 1200 400]);
% term_cmap = [223,194,125; 125,125,125]./255; %shades of brown & gray
term_cmap = [200,200,200; 0,0,0]./255; %shades of gray
%fill the background with seasonal colors
fill([seasons(1,1),13,13,seasons(1,1),seasons(1,1)],[-1,-1,1,1,-1],seas_cmap(1,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
fill([0,seasons(2,1),seasons(2,1),0,0],[-1,-1,1,1,-1],seas_cmap(1,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
fill([seasons(2,1),seasons(3,1),seasons(3,1),seasons(2,1),seasons(2,1)],[-1,-1,1,1,-1],seas_cmap(2,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
fill([seasons(3,1),seasons(4,1),seasons(4,1),seasons(3,1),seasons(3,1)],[-1,-1,1,1,-1],seas_cmap(3,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
fill([seasons(4,1),seasons(1,1),seasons(1,1),seasons(4,1),seasons(4,1)],[-1,-1,1,1,-1],seas_cmap(4,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
%plot the data on top
pl_ref = 1;
for j = 1:length(MP)
    %decide the color for the line based on when it starts to retreat
    if nanmean(MP(j).T.termanom_mo(3:4)) > nanmean(MP(j).T.termanom_mo(5:6))
        cmap_ind = 1; %early retreater!
    else
        cmap_ind = 2; %summer retreater
    end

    % plot based on site size & timing of seasonal retreat
    if ismember(MP(j).name,big3)
        if cmap_ind == 1
            plot([1:12]+0.5,MP(j).T.termanom_mo./max(abs(MP(j).T.termanom_mo)),'-','color',(pl_ref/5)*term_cmap(cmap_ind,:),'linewidth',3); hold on;
            pl(pl_ref) = plot([1:12]+0.5,MP(j).T.termanom_mo./max(abs(MP(j).T.termanom_mo)),'d',...
                'color',(pl_ref/5)*term_cmap(cmap_ind,:),'linewidth',1,'markerfacecolor',(pl_ref/5)*term_cmap(cmap_ind,:)); hold on;
            site_ind(pl_ref) = strmatch(sitenames(j,:),geo_order); pl_ref = pl_ref+1;
        else
            plot([1:12]+0.5,MP(j).T.termanom_mo./max(abs(MP(j).T.termanom_mo)),'--','color',term_cmap(cmap_ind,:),'linewidth',1); hold on;
            plot([1:12]+0.5,MP(j).T.termanom_mo./max(abs(MP(j).T.termanom_mo)),'d',...
                'color',term_cmap(cmap_ind,:),'linewidth',1,'markerfacecolor',term_cmap(cmap_ind,:)); hold on;
        end
    else
        if cmap_ind == 1
            plot([1:12]+0.5,MP(j).T.termanom_mo./max(abs(MP(j).T.termanom_mo)),'-','color',(pl_ref/5)*term_cmap(cmap_ind,:),'linewidth',3); hold on;
            pl(pl_ref) = plot([1:12]+0.5,MP(j).T.termanom_mo./max(abs(MP(j).T.termanom_mo)),'s',...
                'color',(pl_ref/5)*term_cmap(cmap_ind,:),'linewidth',1,'markerfacecolor',(pl_ref/5)*term_cmap(cmap_ind,:)); hold on;
            site_ind(pl_ref) = strmatch(sitenames(j,:),geo_order); pl_ref = pl_ref+1;
        else
            plot([1:12]+0.5,MP(j).T.termanom_mo./max(abs(MP(j).T.termanom_mo)),'--','color',term_cmap(cmap_ind,:),'linewidth',1); hold on;
            plot([1:12]+0.5,MP(j).T.termanom_mo./max(abs(MP(j).T.termanom_mo)),'s',...
                'color',term_cmap(cmap_ind,:),'linewidth',1,'markerfacecolor',term_cmap(cmap_ind,:)); hold on;
        end
    end

end
term_leg = legend(pl,geo_names(site_ind));
set(gca,'xlim',[1,13],'xtick',[1:12],'fontsize',16); grid on;
xlabel('Month','fontsize',16); ylabel('Normalized seasonal terminus anomaly','fontsize',16);
saveas(seasterm_fig,[root_dir,'GrIS-terminus-seasonal-anomalies_plot.png'],'png'); %save the plots
exportgraphics(seasterm_fig,[root_dir,'GrIS-terminus-seasonal-anomalies_plot.tif'],Resolution=600);
close all; 
clear pl term_cmap site_ind;
disp('Done plotting terminus timeseries');


%% extract melange attributes, estimate buttressing, & make overview plots
close all; drawnow;

%set up site-aggregated plots to be arranged somewhat geographically but
%with a hole for the site map to be inserted
% rows = 9; cols = 2; %9 sites in west Greenland
% plot_locs = [1,3,5,7,9,11,13,15,17,18,16,14,12,2];
rows = 7; cols = 3; %9 sites in west Greenland
plot_locs = [2,1,4,7,10,13,16,19,20,21,18,15,12,3];

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
Hfig = figure; set(Hfig,'position',[-1650 50 900 1200]);
% Vfig = figure; set(Vfig,'position',[150 50 1200 1200]);
for j = 1:length(plot_locs)
    figure(Hfig);
    eval(['subH',num2str(geo_ind(j)),'=subplot(',num2str(rows),',',num2str(cols),',',num2str(plot_locs(j)),');']);
    % figure(Vfig);
    % eval(['subV',num2str(geo_ind(j)),'=subplot(',num2str(rows),',',num2str(cols),',',num2str(plot_locs(j)),');']);
end
missfig = figure; set(missfig,'position',[1050 50 450 500]); 
% subm1 = subplot(2,1,1); subm2 = subplot(2,1,2);

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
    %create seasonal average distributions for each fjord bin
    % for k = 1:4
    %     bergdist_seas(k,:,1) = nanmean(bergdist_seaward_setdx_norm(ismember(bmos,seasons(k,:))==1,:),1); %could change to bergdist_seaward_end_norm for moving end
    %     bergdist_seas(k,:,2) = nanmean(bergdist_inland_norm(ismember(bmos,seasons(k,:))==1,:),1);
    % end
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
            
            % %stack the profiles
            % n1_seasmean(:,:,k) = nanmean(binned_n1_stack(:,:,ismember(bmos,seasons(k,:))==1),3);
            % n1_seasmax(:,:,k) = (nanmean(binned_n1_stack(:,:,ismember(bmos,seasons(k,:))==1),3)+std(binned_n1_stack(:,:,ismember(bmos,seasons(k,:))==1),0,3,'omitnan'));
            % n1_seasmin(:,:,k) = (nanmean(binned_n1_stack(:,:,ismember(bmos,seasons(k,:))==1),3)-std(binned_n1_stack(:,:,ismember(bmos,seasons(k,:))==1),0,3,'omitnan'));

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
    % disp([MP(j).name,': j = ',num2str(j),', plot_ind = ',num2str(plot_ind)])

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
            %add symbols
            if ismember(MP(j).name,big3)
                plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'d','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            else
                plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'s','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            end


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
                pH(k) = plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
                if ismember(MP(j).name,big3)
                    % overlay diamond symbols
                    plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'d','color',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'markersize',5); hold on;
                else
                    % overlay square symbols
                    plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'s','color',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'markersize',5); hold on;
                end
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
            'ylim',[0 ceil(max(H_ylim(:,2))/10)*10],'fontsize',12); drawnow;
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
            xlabel('Distance from terminus (km)','fontsize',12); 
            if plot_locs(plot_ind) == rows*cols - (cols-1)
                ax1.YAxis(1).Label.String = 'Thickness (m)';
                ax1.YAxis(1).Label.Position = [-2.7066e3 250 -1];
                % ylabel('Thickness (m)','fontsize',12);
            end
        end
        clear ylims yticks;
        %add a legend
        if j == length(MP)
            H_leg = legend(pH,season_names);
            subH_pos = get(gca,'position');
            H_leg.Orientation = 'horizontal'; %move the legend
            H_leg.Position = [0.38 0.24 0.25 0.035]; 
            set(gca,'position',subH_pos); %resize the subplot
            drawnow;
        end

        %slightly shift plot locations to fit inset map in middle
        if plot_ind > 1 && plot_ind <= 8
            set(gca,'position',[gca_pos(1)-0.06 gca_pos(2) gca_pos(3)+0.04 0.07]);
        elseif plot_ind > 9
            set(gca,'position',[gca_pos(1)+0 gca_pos(2) gca_pos(3)+0.04 0.07]);
        else
            set(gca,'position',[gca_pos(1)-0.03 gca_pos(2) gca_pos(3)+0.04 0.07]);
        end
    end
    
    %add speed profiles to figures
    for k = 1:4
        if k <= 2
            fill_alpha = 0.1;
        else
            fill_alpha = 0.2;
        end

        %relative to moving terminus
        v_mean = nanmean(vel_seas(:,k,:),3)./365;
        vmax = (nanmean(vel_seas(:,k,:),3)+std(vel_seas(:,k,:),0,3,'omitnan'))./365; 
        vmin = (nanmean(vel_seas(:,k,:),3)-std(vel_seas(:,k,:),0,3,'omitnan'))./365; 
        
        %mask out values near zero because they are so much slower than the
        %melange velocities that they must be bad sea ice tracking
        bottomout = find(v_mean < 0.5,1,'first');
        v_mean(bottomout:end) = NaN; vmax(bottomout:end) = NaN; vmin(bottomout:end) = NaN; 

        %plot
        vmax_idx = find(~isnan(vmax)==1); vmin_idx = find(~isnan(vmin)==1);
        if sum(~isnan(v_mean)) ~= 0
            figure(sitefig); subplot(subv);
            fill([vdist(vmax_idx), fliplr(vdist(vmin_idx))]',[vmax(vmax_idx); flipud(vmin(vmin_idx))],...
                seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
            pv(k) = plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
            %add symbols
            if ismember(MP(j).name,big3)
                plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'d','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            else
                plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'s','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            end

            %regional fig
            figure(Hfig);
            % figure(Vfig);
            if ~isempty(plot_ind)
                %navigate to the subplot
                eval(['subplot(subH',num2str(j),');']);
                % eval(['subplot(subV',num2str(j),');']);

                %plot the speed profiles
                yyaxis right;
                fill([vdist(vmax_idx), fliplr(vdist(vmin_idx))]',[vmax(vmax_idx); flipud(vmin(vmin_idx))],...
                    seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
                if k >= 2
                    pV(k-1) = plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
                else
                    plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
                end
                %add symbols
                if ismember(MP(j).name,big3)
                    %overlay diamond symbols
                    plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'d','color',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'markersize',5); hold on;
                else
                    %overlay square symbols
                    plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'s','color',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'markersize',5); hold on;
                end
            end
            MP(j).V.Vseas(k,:) = v_mean';
        else
            figure(sitefig); subplot(subv);
            MP(j).V.Vseas(k,:) = NaN(1,size(vel_seas,1));
        end
        clear v_mean vmax* vmin*;
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
        set(gca,'xlim',[0,25000],'xticklabel',[],'fontsize',12); drawnow;
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
            xlabel('Distance from terminus (km)','fontsize',12);
            if plot_locs(plot_ind) == rows*cols
                ax1.YAxis(2).Label.String = 'Speed (m/d)';
                ax1.YAxis(2).Label.Position = [2.7559e4 550 -1];
                % ylabel('Speed (m/d)','fontsize',12);
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
        loglog(MP(j).D.area,bergdist_seas(k,:,1),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
    end
    set(gca,'ylim',[10^-10 10^-1],'xlim',[10^1 10^7],...
        'ytick',[10^-9,10^-5,10^-1],...
        'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',16); grid on;
    xlabel('Surface area (m^2)','fontsize',16); ylabel('Iceberg count','fontsize',16);
    % text(10000,10^-2,'near-terminus','fontsize',16);
    %seaward bin
    subplot(4,2,6);
    for k = 1:4
        loglog(MP(j).D.area,bergdist_seas(k,:,2),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
    end
    set(gca,'ylim',[10^-10 10^-1],'xlim',[10^1 10^7],...
        'ytick',[10^-9,10^-5,10^-1],...
        'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',16); grid on;
    xlabel('Surface area (m^2)','fontsize',16);
    % text(10000,10^-2,'seaward','fontsize',16);
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
        if sum(~isnan(size_plslope(k,:))) > 0
           ps(k) = plot(Hdist,res_A(k,:),'-','color',seas_cmap(k,:),'linewidth',3); hold on; %bergy bit misfits
           %add symbols
           if ismember(MP(j).name,big3)
               plot(Hdist,res_A(k,:),'d','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on; %bergy bit misfits
           else
               plot(Hdist,res_A(k,:),'s','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on; %bergy bit misfits
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
            %add symbols
           if ismember(MP(j).name,big3)
               plot(Hdist(~isnan(Hmean))',packmean(~isnan(Hmean)),'d','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on; %bergy bit misfits
           else
               plot(Hdist(~isnan(Hmean))',packmean(~isnan(Hmean)),'s','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on; %bergy bit misfits
           end
            ax2.YLim = [0,100]; 
        end

        %format the site figure
        ax2.YAxis(1).Color = 'k'; ax2.YAxis(2).Color = 'k'; 
        ax2.YAxis(1).Label.String = 'Bergy bit misprediction (%)';
        ax2.YAxis(2).Label.String = 'Packing density (%)';

        %scatterplot of bergy bit misfits vs packing density for all sites
        figure(missfig); set(gca,'box','on'); Hdist(Hdist==0) = 1;
        if ismember(MP(j).name,big3)
            %plot dummy points for the legend
            small_ref = find(Hdist<1000,1,'first');
            if ~isempty(small_ref) && ~isnan(packmean(small_ref))
                pm(1) = plot(res_A(k,small_ref),packmean(small_ref),'d','MarkerSize',sqrt(240-18*log(Hdist(small_ref))),'MarkerFaceColor','none',...
                    'MarkerEdgeColor','k','LineWidth',1.5); hold on;
            end
            clear small_ref;
            med_ref = find(Hdist<10000,1,'last');
            if ~isempty(med_ref) && ~isnan(packmean(med_ref))
                pm(2) = plot(res_A(k,med_ref),packmean(med_ref),'d','MarkerSize',sqrt(240-18*log(Hdist(med_ref))),'MarkerFaceColor','none',...
                    'MarkerEdgeColor','k','LineWidth',1.5); hold on;
            end
            clear med_ref;
            big_ref = find(Hdist<20000,1,'last');
            if ~isempty(big_ref)
                pm(3) = plot(res_A(k,big_ref),packmean(big_ref),'d','MarkerSize',sqrt(240-18*log(Hdist(big_ref))),'MarkerFaceColor','none',...
                    'MarkerEdgeColor','k','LineWidth',1.5); hold on;
            end
            clear big_ref;

            %plot all data along the profiles
            scatter(res_A(k,:),packmean,240-18*log(Hdist),'d','MarkerFaceColor','none',...
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
            scatter(res_A(k,:),packmean,240-18*log(Hdist),'s','MarkerFaceColor','none',...
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
        % disp(['Min bergy bit plot symbol size: ',num2str(min(240-18*log(Hdist)))]);
        drawnow; clear packmean Hmean;
    end
    %call out the data for the example map & size distribution plot
    if contains(MP(j).name,'ASG')
        Dsubs = dir([root_dir,MP(j).name,'/',MP(j).name,'*-iceberg-distribution-subsets.csv']);
        for p = 1:length(Dsubs)
            distdate(p,:) = Dsubs(p).name(5:12);
        end
        dateref = strmatch('20110610',distdate);
        %load the data for that date
        D = readtable([root_dir,MP(j).name,'/',Dsubs(dateref).name],"VariableNamingRule","preserve");
        v1 = table2array(D(:,1)); dv1 = table2array(D(:,2));
        %first 2 columns are area & bin width, so 3+ are data from points along the centerline
        if ~isnan(inland_idx(p))
            berg_nos = table2array(D(:,(inland_idx(dateref)-1)+2));
            berg_nos(berg_nos==0) = NaN;
            n1 = berg_nos./dv1;
        end
        clear D;
        clear berg_nos Dsubs distdate;
        packing_dated = packing(dateref,inland_idx(dateref)-1);

        %filter to decide if there are enough points to fit lines
        v = v1(n1>nthresh); n = n1(n1>nthresh); dv = dv1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
        n = n(v>vthresh); dv = dv(v>vthresh); v = v(v>vthresh); %remove smallest size classes
        v = v(~isnan(n)); dv = dv(~isnan(n)); n = n(~isnan(n)); % Remove NaNs

        %enter unfiltered data into model (filtering happens there)
        tempfig = figure;
        [alpha,c1,c2,c3,c4,~,~] = EBC_fragmentation_curve(v1, n1, dv1, norm_type, vthresh, nthresh, 1, normalize_exp); % fit Eq. 1
        n1_mod = EBC_model([c1,c2,alpha,c3,c4],v1); % grab model
        close(tempfig);

        %calculate misfit for bergy bits
        res = n1_mod-n1; % calculate residuals from E-BC model (data - model)
        res_Adated = sum(v1(v1<=vthresh).*(res(v1<=vthresh).*dv1(v1<=vthresh)))./sum(v1.*(n1.*dv1),"omitnan");
        clear v dv n v1 dv1 n1 dateref res;

        %add example date & site info to the scatterplot
        figure(missfig); 
        scatter(res_Adated,100*packing_dated,240,'s','MarkerFaceColor',seas_cmap(2,:),...
            'MarkerEdgeColor',seas_cmap(2,:),'LineWidth',1.5); hold on;
        clear *_dated;
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
    set(gca,'ylim',[-1,1],'ytick',[-1:0.5:1],'yticklabel',[-100:50:100]); %adjust labeling to percents
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
    set(gca,'fontsize',12); grid on; drawnow;
    set(gca,'xlim',[-1,1],'xtick',[-1:0.5:1],'xticklabel',[-100:50:100],'ylim',[0,100]);
    xlabel('Bergy bit misprediction (% Area)','fontsize',12); ylabel('Packing density (%)','fontsize',12); 
    clear pz ps *_leg pos;

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

    %save the data and the figure
    save([root_dir,'GrIS-melange-characteristics.mat'],'MP','-v7.3');
    saveas(sitefig,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-speed-size_',num2str(Hcutoff),'m-Hthreshold_',num2str(vdtmin),'-',num2str(vdtmax),'dt-',vfilter,'-speeds_',sampling,'-profiles.png'],'png'); %save the plots
    exportgraphics(sitefig,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-speed-size_',num2str(Hcutoff),'m-Hthreshold_',num2str(vdtmin),'-',num2str(vdtmax),'dt-',vfilter,'-speeds_',sampling,'-profiles.tif'],Resolution=600);
    % saveas(bergAfig,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-melange-properties_profiles.png'],'png'); %save the plots
    disp(['Done with #',num2str(j),': ',MP(j).name]);
    % % figure(sitefig);
    % % disp('Close the site figure to advance'); disp(' ');
    % % uiwait %advance only after figure is closed
    % saveas(missfig,[root_dir,MP(j).name,'/',MP(j).name,'-bergybit_scatterplot.png'],'png'); 
    % exportgraphics(missfig,[root_dir,MP(j).name,'/',MP(j).name,'-bergybit_scatterplot.tif'],Resolution=600);
    % clf(missfig); %clear everything from the bergy bit misfit figure so you can look at data for each site
    close(sitefig); %close(bergAfig);
    drawnow; 

    %refresh
    clear berg_* bergdist* berg_normdist* C centerline* D Dsubs *idx seaward_* inland_* term_* tran_* size_classes Zfilt H_* Havg vel_* vels* v_mean w zdate berg_mo bins bin_no z_* pos pz pv *dist *yrs *mos seas_leg packing;
    clear leg_* ref_* Tdate* xlims Tslopes res* Tres size_plslope bin* column_names dv1 v1 n1* *_seas ax ax2 max_xlim width_prof Vterm dVdxterm subv subz*;
    clear ax1 ax2;

end
figure(missfig); miss_leg = legend(pm,'terminus','10 km','20 km'); set(miss_leg,'location','southeast');
saveas(missfig,[root_dir,'GrIS-bergybit_scatterplot.png'],'png'); 
exportgraphics(missfig,[root_dir,'GrIS-bergybit_scatterplot.tif'],Resolution=600);
clear pm;


%save the GrIS-wide profiles
saveas(Hfig,[root_dir,'GrIS-melange_thickness-speed_profiles.png'],'png'); 
exportgraphics(Hfig,[root_dir,'GrIS-melange_thickness-speed_profiles.tif'],Resolution=600);
% saveas(Hfig,[root_dir,'GrIS-melange-thickness_profiles.png'],'png'); 
% saveas(Vfig,[root_dir,'GrIS-melange-speed_profiles.png'],'png'); 


%create histograms of buttressing
buttfig = figure; set(buttfig,'position',[450 50 900 900]); 
for k = [2,3,4,1]
    subplot(2,2,1);
    hM(k) = histogram(BM_annual(~isnan(BM_annual(:,k)),k),'BinEdges',[0.1:0.1:12],'FaceColor',seas_cmap(k,:),'EdgeColor',seas_cmap(k,:),...
        'EdgeAlpha',1,'LineWidth',1); hold on;
    set(gca,'fontsize',16,'ylim',[0,17],'xlim',[0,12],'box','on'); ylabel('Count','fontsize',16); xlabel('Buttressing (x10^6 N/m)','fontsize',16);
    text(0.15,0.95*17,'a) annual packing-based buttressing','fontsize',16);
    subplot(2,2,2);
    hA(k) = histogram(BA_annual(~isnan(BA_annual(:,k)),k),'BinEdges',[0.1:0.1:12],'FaceColor',seas_cmap(k,:),'EdgeColor',seas_cmap(k,:),...
        'EdgeAlpha',1,'LineWidth',1); hold on;
    set(gca,'fontsize',16,'ylim',[0,17],'xlim',[0,12],'box','on'); xlabel('Buttressing (x10^6 N/m)','fontsize',16);
    text(0.15,0.95*17,'b) annual strain rate-based buttressing','fontsize',16);
    subplot(2,2,3);
    histogram(BM_character(~isnan(BM_character(:,k)),k),'BinEdges',[0.1:0.1:7],'FaceColor',seas_cmap(k,:),'EdgeColor',seas_cmap(k,:),...
        'EdgeAlpha',1,'LineWidth',1); hold on;
    set(gca,'fontsize',16,'ylim',[0,6],'xlim',[0,7],'box','on'); ylabel('Count','fontsize',16); xlabel('Buttressing (x10^6 N/m)','fontsize',16);
    text(0.1,0.95*6,'c) characteristic packing-based buttressing','fontsize',16);
    subplot(2,2,4);
    histogram(BA_character(~isnan(BA_character(:,k)),k),'BinEdges',[0.1:0.1:7],'FaceColor',seas_cmap(k,:),'EdgeColor',seas_cmap(k,:),...
        'EdgeAlpha',1,'LineWidth',1); hold on;
    set(gca,'fontsize',16,'ylim',[0,6],'xlim',[0,7],'box','on'); xlabel('Buttressing (x10^6 N/m)','fontsize',16);
    text(0.1,0.95*6,'d) characteristic strain rate-based buttressing','fontsize',16);

    %display medians across all sites for the season
    disp(['Median Meng buttressing estimate for ',char(season_names(k)),' = ',num2str(round(nanmedian(BM_character(:,k),"all"),2)),'x10^6 N/m']);
    disp(['Median Amundson buttressing estimate for ',char(season_names(k)),' = ',num2str(round(nanmedian(BA_character(:,k),"all"),2)),'x10^6 N/m']);
end
subplot(2,2,1); pos = get(gca,'position'); set(gca,'position',[pos(1)-0.03 pos(2) 1.15*pos(3) 1.05*pos(4)]);
subplot(2,2,2); pos = get(gca,'position'); set(gca,'position',[pos(1)-0.03 pos(2) 1.15*pos(3) 1.05*pos(4)]);
subplot(2,2,3); pos = get(gca,'position'); set(gca,'position',[pos(1)-0.03 pos(2) 1.15*pos(3) 1.05*pos(4)]);
subplot(2,2,4); pos = get(gca,'position'); set(gca,'position',[pos(1)-0.03 pos(2) 1.15*pos(3) 1.05*pos(4)]);
disp(['Mean difference between Meng & Amundson buttressing estimates = ',num2str(round(nanmean(BM_annual-BA_annual,"all"),2)),'x10^6 N/m']);
disp(['Mean PERCENT difference between Meng & Amundson buttressing estimates = ',num2str(round(100*nanmean((BM_annual-BA_annual)./nanmean(cat(3,BM_annual,BA_annual),3),"all"),0)),'%']);
%display medians across spring & summer using only sites with data from both seasons
disp('Seasonal Meng & Amundson buttressing omitting sites with spring-only data:')
disp(['Medians for ',char(season_names(2)),' = ',num2str(round(nanmedian(BM_character(~isnan(BM_character(:,3)),2),"all"),2)),...
    ' and ',num2str(round(nanmedian(BA_character(:,2),"all"),2)),'x10^6 N/m']);
disp(['Medians for ',char(season_names(3)),' = ',num2str(round(nanmedian(BM_character(~isnan(BM_character(:,3)),3),"all"),2)),...
    ' and ',num2str(round(nanmedian([BA_character(~isnan(BA_character(:,3)),3); zeros(sum(~isnan(BM_character(~isnan(BM_character(:,3)),2)))-sum(~isnan(BA_character(:,3))),1)],"all"),2)),'x10^6 N/m']);
disp('Seasonal Meng & Amundson buttressing omitting sites with spring-only data AND without strain rate-based estimates:')
disp(['Medians for ',char(season_names(2)),' = ',num2str(round(nanmedian(BM_character(~isnan(BM_character(:,3)) & ~isnan(BA_character(:,3)),2),"all"),2)),...
    ' and ',num2str(round(nanmedian(BA_character(~isnan(BM_character(:,3)) & ~isnan(BA_character(:,3)),2),"all"),2)),'x10^6 N/m']);
disp(['Medians for ',char(season_names(3)),' = ',num2str(round(nanmedian(BM_character(~isnan(BM_character(:,3)) & ~isnan(BA_character(:,3)),3),"all"),2)),...
    ' and ',num2str(round(nanmedian(BA_character(~isnan(BM_character(:,3)) & ~isnan(BA_character(:,3)),3),"all"),2)),'x10^6 N/m']);
%save the buttressing plot
saveas(buttfig,[root_dir,'GrIS-melange_buttressing_histograms.png'],'png'); 
exportgraphics(buttfig,[root_dir,'GrIS-melange_buttressing_histograms.tif'],Resolution=600);

%export overview data to a CSV (used for plotting sitemap in QGIS)
TS = [array2table(site_naming(:,2)), array2table(PS_X), array2table(PS_Y), array2table(sp_stats)];
column_names = ["Site Name","X (m)", "Y (m)","Melange thickness (m)","Melange buttressing (x10^6 N/m)"]; TS.Properties.VariableNames = column_names;
writetable(TS,[root_dir,'GrIS-melange-sites.csv']);
clear TS;

%% display the seasonal information for each site (GRL paper Table)
disp('Seasonal statistics for all sites:')

for j = 1:length(geo_ind)
    %grab date info
    for p = 1:length(MP(geo_ind(j)).Z.date)
        zdate(p) = convert_to_decimaldate(char(MP(geo_ind(j)).Z.date(p)));
        datest(p,:) = datetime(MP(geo_ind(j)).Z.date{p},'InputFormat','yyyyMMdd');
        dateout(p,:) = datestr(datest(p,:),'yyyy-mm-dd');
        zyrs(p) = year(datest(p,:)); zmos(p) = month(datest(p,:));
        % clear datest;
    end

    %site abbreviation and name
    disp(MP(geo_ind(j)).name);
    disp(char(geo_names(j)));

    %seasonal DEM dates, near-terminus thickness, buttressing
    for k = 1:4
        disp(['  ',char(season_names(k))]);
        %DEM dates
        seas_refs = find(ismember(zmos,seasons(k,:))==1);
        if ~isempty(seas_refs)
            date_cat = ['    '];
            for p = 1:length(seas_refs)
                if p ~= 1
                    date_cat = [date_cat,', ',dateout(seas_refs(p),:)];
                else
                    date_cat = [date_cat,dateout(seas_refs(p),:)];
                end
            end
            disp(date_cat)
        end
        clear seas_refs;

        %near-terminus thickness
        disp(['    thickness (m): ',num2str(round(nanmean(MP(geo_ind(j)).B.Ho(zcutoff+1,k,:)),1))]);

        %buttressing
        disp(['    packing-based buttressing (N/m): ',num2str(round(nanmean(MP(geo_ind(j)).B.butt_Meng(zcutoff+1,k,:))./10^6,2))]);
        disp(['    strainrate-based buttressing (N/m): ',num2str(round(nanmean(MP(geo_ind(j)).B.butt_Amundson(zcutoff+1,k,:))./10^6,2))]);
    end

    clear zdate datest dateout zyrs zmos;
end


%% create overview maps for each site
cd(root_dir);
load([root_dir,'GrIS-melange-characteristics.mat']);

%plot option
suffix = 'site-map_bigfont'; %add 'bigfont' to the suffix if using 30pt
map_font = 28; %default font = 16, bigfont = 28
if contains(suffix,'bigfont')
    x_shift = 0.05;
else
    x_shift = 0.02;
end

%iterate
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
    %crop the image to adjust brightnesses appropriately
    xlims = [find(im.x<=min(melmask.uncropped.x),1,'last')-500/R.CellExtentInWorldX, find(im.x<=max(melmask.uncropped.x),1,'last')+500/R.CellExtentInWorldX];
    ylims = [find(im.y>=max(melmask.uncropped.y),1,'last')-500/R.CellExtentInWorldY, find(im.y<=min(melmask.uncropped.y),1,'first')+500/R.CellExtentInWorldY];
    xlims(xlims<1) = 1; xlims(xlims>length(im.x)) = length(im.x);
    ylims(ylims<1) = 1; ylims(ylims>length(im.y)) = length(im.y);
    im_subset = im.z(min(ylims):max(ylims),min(xlims):max(xlims),:);
    im_subset = im_subset./max(max(im_subset));
    clear I R;

    %isolate terminus positions from the dated melange masks
    for p = 1:length(melmask.dated)
        [~,on] = inpolygon(melmask.dated(p).x,melmask.dated(p).y,melmask.uncropped.x,melmask.uncropped.y);
        mask_xrange(p,:) = [min(melmask.dated(p).x),max(melmask.dated(p).x)]; 
        mask_yrange(p,:) = [min(melmask.dated(p).y),max(melmask.dated(p).y)];
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
    map_fig = figure; set(map_fig,'position',[850 50 900 800]); ax1 = axes;
    % imagesc(ax1,im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset,[],[])); axis xy equal; 
    imagesc(ax1,im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),im_subset); axis xy equal; 
    colormap(ax1,'gray'); drawnow; hold on;
    set(ax1,'xlim',[min(mask_xrange,[],'all')-500 max(mask_xrange,[],'all')+500],'ylim',[min(mask_yrange,[],'all')-500 max(mask_yrange,[],'all')+500]);
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
    % elev_cmap(1,:) = [1 1 1]; %white

    %plot the time-averaged DEM & redraw the melange mask
    % figure1 = figure; set(gcf,'position',[50 50 1600 600]);
    ax2 = axes; 
    implot = imagesc(ax2,imx,imy,melange_obs/p,'AlphaData', meltrans); axis xy equal; hold on;
    colormap(ax2,elev_cmap); set(ax2,'color','none','visible','off');
    cbar = colorbar; cbar.YLabel.String = 'Iceberg cover frequency';
    cbar.FontSize = map_font-4; cbar.Location = "eastoutside";
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
    plot(ax2,melmask.dated(minind).x,melmask.dated(minind).y,'-','color',[90,174,97]/255,'linewidth',3); hold on;
    plot(ax2,melmask.dated(maxind).x,melmask.dated(maxind).y,'-','color',[90,174,97]/255,'linewidth',3); hold on;
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
    plot(ax2,MP(j).V.X(min(seaward_ext):max(inland_ext)+1),MP(j).V.Y(min(seaward_ext):max(inland_ext)+1),'-k','linewidth',2); hold on;
    if contains(big3,MP(j).name)
        plot(ax2,MP(j).V.X(min(seaward_ext):max(inland_ext)+1),MP(j).V.Y(min(seaward_ext):max(inland_ext)+1),'dk','markersize',6,'linewidth',1,'markerfacecolor','k'); hold on;
    else
        plot(ax2,MP(j).V.X(min(seaward_ext):max(inland_ext)+1),MP(j).V.Y(min(seaward_ext):max(inland_ext)+1),'sk','markersize',6,'linewidth',1,'markerfacecolor','k'); hold on;
    end
    set(ax2,'xlim',[min(mask_xrange,[],'all')-500 max(mask_xrange,[],'all')+500],'ylim',[min(mask_yrange,[],'all')-500 max(mask_yrange,[],'all')+500]);
    % xlims = [min(MP(j).Z.transectX(min(seaward_ext):max(inland_ext)+1,:),[],'all'), max(MP(j).Z.transectX(min(seaward_ext):max(inland_ext)+1,:),[],'all')];
    % ylims = [min(MP(j).Z.transectY(min(seaward_ext):max(inland_ext)+1,:),[],'all'), max(MP(j).Z.transectY(min(seaward_ext):max(inland_ext)+1,:),[],'all')];
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    % set(ax1,'xlim',xlims,'ylim',ylims);
    xticks = get(ax1,'xtick'); yticks = get(ax1,'ytick');
    set(ax1,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',map_font);
    xlabel(ax1,'Easting (km)','fontsize',map_font); ylabel(ax1,'Northing (km)','fontsize',map_font);
    pos = get(gca,'position');
    if range(xlims) > 1.05*range(ylims) %short and fat map so plot the legend below
        % map_leg = legend(pt,season_names,'Location','southoutside',...
        %     'Orientation','horizontal');
        % leg_pos = get(map_leg,'position');
        % set(map_leg,'position',[leg_pos(1)+(0.5-mean([leg_pos(1) leg_pos(1)+leg_pos(3)])) 0.075 leg_pos(3) leg_pos(4)]);
        cbar.Location = "northoutside"; cbar.Orientation = 'horizontal';

        % if range(xlims)/range(ylims) < 1.5
            set(ax1,'position',[pos(1)+x_shift 0.2 0.9*pos(3) 0.9*pos(4)]); drawnow;
            set(ax2,'position',[pos(1)+x_shift 0.2 0.9*pos(3) 0.9*pos(4)]); drawnow;
        % else
        %     set(ax1,'position',[pos(1)+x_shift 0.2 pos(3) pos(4)]); drawnow;
        %     set(ax2,'position',[pos(1)+x_shift 0.2 pos(3) pos(4)]); drawnow;
        % end
    else %tall and thin map so plot the legend on the side
        % map_leg = legend(pt,season_names,'Location','eastoutside',...
        %     'Orientation','vertical');
        % leg_pos = get(map_leg,'position');
        % if leg_pos(1)+0.05 > pos(1)-0.05+pos(3)
        %     set(map_leg,'position',[leg_pos(1)+0.05 leg_pos(2) leg_pos(3) leg_pos(4)]);
        % else
        %     set(map_leg,'position',[pos(1)-0.05+pos(3) leg_pos(2) leg_pos(3) leg_pos(4)]);
        % end
        set(ax1,'position',[pos(1) pos(2) pos(3) pos(4)]); drawnow;
        set(ax2,'position',[pos(1) pos(2) pos(3) pos(4)]); drawnow;
    end
    saveas(map_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-',suffix,'.png'],'png'); %save the image
    exportgraphics(map_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-',suffix,'.tif'],Resolution=600);
    % uiwait %advance only after figure is closed

    %for Alison Glacier for 20110610, create a two subpanel figure to
    %demonstrate how packing density and the bergy bit misfit are
    %calculated with the map on the top & the size distribution on the bottom
    if contains(MP(j).name,'ASG')
        ex_fig = figure; set(ex_fig,'position',[950 50 450 500]); 
        subm = subplot(2,1,1); 

        %load the DEM of interest
        melange_mats = dir([root_dir,site_abbrev,'/DEMs/*DEMfilled.mat']);
        for p = 1:length(melange_mats)
            DEM_name = melange_mats(p).name;
            if contains(DEM_name,'20110610')
                load([root_dir,site_abbrev,'/DEMs/',DEM_name]);
                [ZXgrid,ZYgrid] = meshgrid(M.DEM.x,M.DEM.y);
                
                %create the packing density mask 
                packmask = zeros(size(ZXgrid));
                maskedz = M.DEM.z_filled.*M.mask.DEM;
                packmask(maskedz <= zthresh) = 1;

                clear Z*;
            end
        end

        %plot the map
        maskedz(isnan(maskedz)) = -0.1; axa = gca;
        imagesc(subm,M.DEM.x,M.DEM.y,maskedz); axis xy equal; 
        berg_cmap = colormap(gray(501)); berg_cmap(1,:) = [1 1 1];
        colormap(subm,berg_cmap); drawnow; hold on;
        % set(subm,'xlim',[min(mask_xrange,[],'all')-500 max(mask_xrange,[],'all')+500],'ylim',[min(mask_yrange,[],'all')-500 max(mask_yrange,[],'all')+500]);
        set(subm,'xlim',[-332500 -320500],'ylim',[-1646500 -1638500]);
        set(subm,'clim',[0 50]); berg_cbar = colorbar; berg_cbar.Ticks = [0:10:50]; berg_cbar.Label.String = 'elevation (m a.s.l.)';
        %format the map
        xticks = get(subm,'xtick'); yticks = get(subm,'ytick');
        set(subm,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',12);
        xlabel(subm,'Easting (km)','fontsize',12); ylabel(subm,'Northing (km)','fontsize',12);
        

        %load the size distribution for near the terminus for the date
        D = readtable([root_dir,MP(j).name,'/',MP(j).name,'-20110610-iceberg-distribution-subsets.csv'],"VariableNamingRule","preserve");
        v1 = table2array(D(:,1)); dv1 = table2array(D(:,2));
        berg_nos = table2array(D(:,3:end)); berg_nos(berg_nos==0) = NaN;
        size_classes = sum(~isnan(berg_nos),1);
        inland_idx = find(size_classes>0,1,'last');
        n1 = berg_nos(:,inland_idx)./dv1;
        clear D;
        %filter the data & fit the curve
        v = v1(n1>nthresh); n = n1(n1>nthresh); dv = dv1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
        n = n(v>vthresh); dv = dv(v>vthresh); v = v(v>vthresh); %remove smallest size classes
        v = v(~isnan(n)); dv = dv(~isnan(n)); n = n(~isnan(n)); % Remove NaNs
        %enter unfiltered data into model (filtering happens there)
        tempfig = figure; drawnow;
        figure(tempfig);
        [alpha,c1,c2,c3,c4,~,~] = EBC_fragmentation_curve(v1, n1, dv1, norm_type, vthresh, nthresh, 1, normalize_exp); % fit Eq. 1
        n1_mod = EBC_model([c1,c2,alpha,c3,c4],v1); % grab model
        close(tempfig);

        %add the mask for bergy bits
        packmask(maskedz==-0.1) = -1; packmask(maskedz==0) = -1;
        packalpha = 0.75*packmask;
        figure(ex_fig); subplot(subm);
        axb = axes; mask_cmap = [1,1,1; 0,0,0; 251/255,180/255,185/255];
        implot = imagesc(axb,M.DEM.x,M.DEM.y,packmask,'AlphaData', packalpha); axis xy equal; hold on;
        colormap(axb,mask_cmap); set(axb,'color','none','visible','off');
        % set(axb,'xlim',[min(mask_xrange,[],'all')-500 max(mask_xrange,[],'all')+500],'ylim',[min(mask_yrange,[],'all')-500 max(mask_yrange,[],'all')+500]);
        set(axb,'xlim',[-332500 -320500],'ylim',[-1646500 -1638500]);
        colormap(subm,berg_cmap);
        map_pos = get(subm,'position'); set(axb,'position',map_pos);
        clear packmask maskedz M;
        % %add the centerline
        % plot(axb,MP(j).V.X(min(seaward_ext):max(inland_ext)+1),MP(j).V.Y(min(seaward_ext):max(inland_ext)+1),'-k','markersize',6,'linewidth',2); hold on;
        % plot(axb,MP(j).V.X(min(seaward_ext):max(inland_ext)+1),MP(j).V.Y(min(seaward_ext):max(inland_ext)+1),'sk','markersize',6,'linewidth',1,'markerfacecolor','k'); hold on;
        %load the transects
        S = shaperead([root_dir,MP(j).name,'/shapefiles/',MP(j).name,'-transects_1000m.shp']);
        clear *_ext;

        %plot the transect bounds for the subset
        % plot(melmask.uncropped.x,melmask.uncropped.y,'-k'); hold on;
        end_dist = sqrt((melmask.uncropped.x-S(inland_idx).X(1)).^2 + (melmask.uncropped.y-S(inland_idx).Y(1)).^2);
        [~,idx] = min(end_dist); end_idx(1) = idx; clear end_dist idx;
        end_dist = sqrt((melmask.uncropped.x-S(inland_idx).X(end-1)).^2 + (melmask.uncropped.y-S(inland_idx).Y(end-1)).^2);
        [~,idx] = min(end_dist); end_idx(2) = idx; clear end_dist idx;
        end_dist = sqrt((melmask.uncropped.x-S(inland_idx+1).X(1)).^2 + (melmask.uncropped.y-S(inland_idx+1).Y(1)).^2);
        [~,idx] = min(end_dist); start_idx(1) = idx; clear end_dist idx;
        end_dist = sqrt((melmask.uncropped.x-S(inland_idx+1).X(end-1)).^2 + (melmask.uncropped.y-S(inland_idx+1).Y(end-1)).^2);
        [~,idx] = min(end_dist); start_idx(2) = idx; clear end_dist idx;
        plot([S(inland_idx).X(1:end-1),melmask.uncropped.x(end_idx(2)+1:1:start_idx(2))',S(inland_idx+1).X(end-1:-1:1),melmask.uncropped.x(start_idx(1)+1:1:end_idx(1)-1)',S(inland_idx).X(1)],...
            [S(inland_idx).Y(1:end-1),melmask.uncropped.y(end_idx(2)+1:1:start_idx(2))',S(inland_idx+1).Y(end-1:-1:1),melmask.uncropped.y(start_idx(1)+1:1:end_idx(1)-1)',S(inland_idx).Y(1)],'-','linewidth',2,...
            'color',seas_cmap(2,:)); hold on;
        clear S *_idx;

        %plot the data
        figure(ex_fig);
        subp = subplot(2,1,2);
        plot(subp,log10([v1(zthresh+1),v1(zthresh+1)]),log10([10^-12 10^3]),'-.','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
        plot(subp,log10(v1),log10(n1),'-','LineWidth',2,'Color',seas_cmap(2,:)); hold on;
        plot(subp,log10(v1),log10(n1),'s','LineWidth',1,'Color',seas_cmap(2,:),'markersize',5,'markerfacecolor',seas_cmap(2,:)); hold on;
        plot(subp,log10(v1),log10(n1_mod),'--','LineWidth',2,'Color',seas_cmap(2,:)); hold on;
        %format the plot
        subp.FontSize = 12;
        subp.YLim = [-5 2]; subp.YTick = [-4:2:2]; %axl.YLim = [10^-10 10^3];
        subp.XLim = [1.5 5]; subp.XTick = [2, 3, 4, 5];
        subp.XTickLabel = compose('%.0e',[10^2, 10^3, 10^4, 10^5]);
        subp.YTickLabel = compose('%.0e', [10^-4, 10^-2, 0, 10^2]);
        xlabel('Surface area (m^2)','fontsize',12);
        ylabel('Normalized iceberg count','fontsize',12);
        plot_pos = get(gca,'position');
        % set(subp,'position',[map_pos(1) plot_pos(2) map_pos(3) map_pos(4)]);
        set(subp,'position',[0.185 0.11 0.7 0.34]);

        %save the figure
        saveas(ex_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-bergy-bit_example-subplots.png'],'png'); %save the image
        exportgraphics(ex_fig,[root_dir,sitenames(j,:),'/',sitenames(j,:),'-bergy-bit_example-subplots.tif'],Resolution=600);
        close(ex_fig);
        clear inland_idx v n dv n1 v1 dv1;
    end

    clear xlims ylims pt elev_cmap maskedz meltrans mel_data termdists maxind minind mask_*range;
    clear berg_* im* LCdir *_ext leg* map_fig map_leg melmask zmos on pos ref* seaward_ext size_classes sort_ind T T_inds term_* Tyrs Tmos *lims zyrs zdate zseas Tseas tran_cmap colorgrad term;
end
close all;
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

%% plot seasonal size distributions with empirical fits to demonstrate bergy bit misfits

%specify if you want to zoom in the plots to focus on the power-law parts
%of the size distributions
suffix = 'plzoom_loglog'; %'plzoom_loglog' to focus on iceberg areas <10^5m^2, 'fulldist_loglog' for full

% %Size distribution parameters (only needed if you didn't run the sectionabove with the same info)
% nthresh = 1e-6; % set small number bin cutoff (n1 must be greater than this value)
% zthresh = 3; %set small size bin cutoff (freeboard must exceed this value)
% rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)
% ARcomp.best.autoALL = 2; % iceberg aspect ratio 
% vthresh = (1/4)*pi*((rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL.*zthresh).^2; %don't include this bin size or smaller in curve fitting
% dplawthresh = 10^5; % upper bound on the intercept for the dummy powerlaw
% norm_type = 2; % toggle between L2, max, and log norm using 2, Inf, and 'log'
% normalize_exp = 1.5; % Increase to weight residuals towards end of the curve, minimum = 1

%create a seasonal profile figure with a top row containing profiles near the terminus and a bottom row containing seaward profiles
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
        bmos(p) = month(berg_dates(p)); byrs(p) = year(berg_dates(p));

        %identify observational limits along the centerline
        berg_nos = table2array(D(:,3:end)); berg_nos(berg_nos==0) = NaN;
        size_classes(p,:) = sum(~isnan(berg_nos),1);
        % seaward_ext(p) = find(size_classes(p,:)>0,1,'first');
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
    clear size_classes;

    %now analyze the size distributions
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

            size_classes = sum(~isnan(berg_nos),1);
            seaward_ext(p) = find(size_classes>0,1,'last'); %seaward limit index in moving reference frame relative to the terminus
        end
        clear D;
        %stack distributions in a structure
        binned_n1_stack(:,:,p) = n1;
        clear n1 berg_nos size_classes; 
    end
    %identify seaward sample location
    seaward_meanidx = round(mean(seaward_ext));


    %pasted from above...
    %create seasonal average distributions for each fjord bin
    n1wi_stack = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
    n1sp_stack = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
    n1su_stack = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
    n1fa_stack = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
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
                n1wi_stack(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(1,:))==1),3);
            end
            if ~isempty(ismember(mos_yr,seasons(2,:))) %spring
                n1sp_stack(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(2,:))==1),3);
            end
            if ~isempty(ismember(mos_yr,seasons(3,:))) %summer
                n1su_stack(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(3,:))==1),3);
            end
            if ~isempty(ismember(mos_yr,seasons(4,:))) %fall
                n1fa_stack(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(4,:))==1),3);
            end

            clear mos_yr n1_temp;
        end
        clear yr_idx
    end
    %calculate seasonal averages across all years & fit curves
    n1_seasmean(:,:,1) = nanmean(n1wi_stack,3);
    n1_seasmax(:,:,1) = (nanmean(n1wi_stack,3)+std(n1wi_stack,0,3,'omitnan')); n1_seasmin(:,:,1) = (nanmean(n1wi_stack,3)-std(n1wi_stack,0,3,'omitnan'));
    n1_seasmean(:,:,2) = nanmean(n1sp_stack,3);
    n1_seasmax(:,:,2) = (nanmean(n1sp_stack,3)+std(n1sp_stack,0,3,'omitnan')); n1_seasmin(:,:,2) = (nanmean(n1sp_stack,3)-std(n1sp_stack,0,3,'omitnan'));
    n1_seasmean(:,:,3) = nanmean(n1su_stack,3);
    n1_seasmax(:,:,3) = (nanmean(n1su_stack,3)+std(n1su_stack,0,3,'omitnan')); n1_seasmin(:,:,3) = (nanmean(n1su_stack,3)-std(n1su_stack,0,3,'omitnan'));
    n1_seasmean(:,:,4) = nanmean(n1fa_stack,3);
    n1_seasmax(:,:,4) = (nanmean(n1fa_stack,3)+std(n1fa_stack,0,3,'omitnan')); n1_seasmin(:,:,4) = (nanmean(n1fa_stack,3)-std(n1fa_stack,0,3,'omitnan'));
    clear n1wi_stack n1sp_stack n1su_stack n1fa_stack;

    %create seasonal average distributions for each fjord bin
    for k = 1:4
        for l = [1,seaward_meanidx]
            n1 = n1_seasmean(:,l,k);
            if l == 1
                n1_inlandseas(:,j,k) = n1;
            else
                n1_seawardseas(:,j,k) = n1;
            end
            clear n1;
        end
    end
    %concatenate the profiles for each season so that they can be plotted
    %as filled polygons instead of individual lines
    %seaward
    n1wi(:,j,2) = n1_seawardseas(:,j,1); n1sp(:,j,2) = n1_seawardseas(:,j,2); 
    n1su(:,j,2) = n1_seawardseas(:,j,3); n1fa(:,j,2) = n1_seawardseas(:,j,4);
    %inland
    n1wi(:,j,1) = n1_inlandseas(:,j,1); n1sp(:,j,1) = n1_inlandseas(:,j,2); 
    n1su(:,j,1) = n1_inlandseas(:,j,3); n1fa(:,j,1) = n1_inlandseas(:,j,4);

    clear berg_dates bmos byrs berg_nos size_classes seaward_* inland_* binned_n1_stack n1_seasmean term_trace;
    clear n1_seasm*;
end
log_v = log10(v1); 

%apply the EBC fitc function to the average of all seasonal profiles in each
%location & plot the result as a dashed line
sample_locs = {'near-terminus','seaward'};
for l = 1:2
    if l == 1
        disp('plotting profiles for near-terminus sample location')
    else
        disp('plotting profiles for seaward sample location')
    end

    %2 = spring, 3 = summer, 4 = fall
    for k = 2:4
        if l == 1
            n1 = nanmedian(n1_inlandseas(:,:,k),2);
        else
            n1 = nanmedian(n1_seawardseas(:,:,k),2);
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
                if k == 2
                    subplot(subia); log_n = log10(n1sp); 
                elseif k == 3
                    subplot(subib); log_n = log10(n1su);
                else
                    subplot(subic); log_n = log10(n1fa);
                end
                yyaxis left; axl = gca; pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);

                %plot dummy fill polygons for the legend
                if l == 1 && k == 2
                    pf(1) = fill(gca,[10,50,50,10,10],[0.9e2,0.9e2,1.1e2,1.1e2,0.9e2],...
                        seas_cmap(2,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
                    pf(2) = fill(gca,[10,50,50,10,10],[0.9e2,0.9e2,1.1e2,1.1e2,0.9e2],...
                        seas_cmap(3,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
                    pf(3) = fill(gca,[10,50,50,10,10],[0.9e2,0.9e2,1.1e2,1.1e2,0.9e2],...
                        seas_cmap(4,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
                end

                %plot the data
                plot(log10([v1(zthresh+1),v1(zthresh+1)]),log10([10^-12 10^3]),'-.','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                plot(log10(v1),log10(n1),'.-','LineWidth',2,'Color','k'); hold on;
                plot(log10(v1),log10(n1_mod),'--','LineWidth',2,'Color','k'); hold on;
                % loglog([v1(zthresh+1),v1(zthresh+1)],[10^-12 10^3],'-','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                % loglog(v1,n1,'.-','LineWidth',2,'Color','k'); hold on;
                % loglog(v1,n1_mod,'--','LineWidth',2,'Color','k'); hold on;

                %format the subplot
                axl.YAxis(1).Color = 'k'; set(axl,'box','on'); axl.LineWidth = 1.5;
                if contains(suffix,'zoom')
                    axl.YLim = [-6 3]; axl.YTick = [-6:2:2]; %axl.YLim = [10^-10 10^3];
                    axl.XLim = [1.5 5]; axl.XTick = [2, 3, 4, 5];
                    axl.XTickLabel = compose('%.0e',[10^2, 10^3, 10^4, 10^5]);
                    if k == 2
                        axl.YTickLabel = compose('%.0e', [10^-6, 10^-4, 10^-2, 0, 10^2]);
                    else
                        axl.YTickLabel = [];
                    end
                    text(3.5,2.3,[alphabet((l-1)*3 + (k-1)),') ',char(sample_locs(l))],'fontsize',16);
                else
                    axl.YLim = [-8 3]; axl.YTick = [-8:2:2]; %axl.YLim = [10^-10 10^3];
                    axl.XLim = [1.5 6.5]; axl.XTick = [2, 4, 6];
                    axl.XTickLabel = compose('%.0e',[10^2, 10^4, 10^6]);
                    if k == 2
                        axl.YTickLabel = compose('%.0e', [10^-8, 10^-6, 10^-4, 10^-2, 0, 10^2]);
                    else
                        axl.YTickLabel = [];
                    end
                    text(4.2,2.2,[alphabet((l-1)*3 + (k-1)),') ',char(sample_locs(l))],'fontsize',16);
                end
                axl.FontSize(1) = 16; %grid on;
                %axl.XLim = [10^1 10^7]; axl.XTick = [10^2, 10^4, 10^6];
                % set(gca,'ylim',[10^-10 10^3],'xlim',[10^1 10^7],'xtick',[10^2, 10^4, 10^6],'fontsize',16); grid on;
                % text(10^5,10^2,[alphabet(k*l-1),') ',char(season_names(k))],'fontsize',16)
                if k == 2
                    ylabel('Normalized iceberg count','fontsize',16);
                end
                
            else
                if k == 2
                    subplot(subsa); log_n = log10(n1sp); 
                elseif k == 3
                    subplot(subsb); log_n = log10(n1su);
                else
                    subplot(subsc); log_n = log10(n1fa);
                end
                yyaxis left; axl = gca; pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);

                %plot the data
                plot(log10([v1(zthresh+1),v1(zthresh+1)]),log10([10^-12 10^3]),'-.','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                plot(log10(v1),log10(n1),'.-','LineWidth',2,'Color','k'); hold on;
                plot(log10(v1),log10(n1_mod),'--','LineWidth',2,'Color','k'); hold on;
                % loglog([v1(zthresh+1),v1(zthresh+1)],[10^-12 10^3],'-','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                % loglog(v1,n1,'.-','LineWidth',2,'Color','k'); hold on;
                % loglog(v1,n1_mod,'--','LineWidth',2,'Color','k'); hold on;

                %format the subplot
                axl.YAxis(1).Color = 'k'; set(axl,'box','on'); axl.LineWidth = 1.5;
                if contains(suffix,'zoom')
                    axl.YLim = [-6 3]; axl.YTick = [-6:2:2]; %axl.YLim = [10^-10 10^3];
                    axl.XLim = [1.5 5]; axl.XTick = [2, 3, 4, 5];
                    axl.XTickLabel = compose('%.0e',[10^2, 10^3, 10^4, 10^5]);
                    if k == 2
                        axl.YTickLabel = compose('%.0e', [10^-6, 10^-4, 10^-2, 0, 10^2]);
                    else
                        axl.YTickLabel = [];
                    end
                    text(3.5,2.3,[alphabet((l-1)*3 + (k-1)),') ',char(sample_locs(l))],'fontsize',16);
                else
                    axl.YAxis(1).Color = 'k'; set(axl,'box','on'); axl.LineWidth = 1.5;
                    axl.YLim = [-8 3]; axl.YTick = [-8:2:2]; %axl.YLim = [10^-10 10^3];
                    axl.XLim = [1.5 6.5]; axl.XTick = [2, 4, 6]; %axl.XLim = [10^1 10^7]; axl.XTick = [10^2, 10^4, 10^6];
                    axl.XTickLabel = compose('%.0e',[10^2, 10^4, 10^6]);
                    if k == 2
                        axl.YTickLabel = compose('%.0e', [10^-8, 10^-6, 10^-4, 10^-2, 0, 10^2]);
                    else
                        axl.YTickLabel = [];
                    end
                    text(4.2,2.2,[alphabet((l-1)*3 + (k-1)),') ',char(sample_locs(l))],'fontsize',16);
                end
                axl.FontSize(1) = 16; %grid on;
                % set(gca,'ylim',[10^-10 10^3],'xlim',[10^1 10^7],'xtick',[10^2, 10^4, 10^6],'fontsize',16); 
                % text(10^5,10^2,[alphabet((k-1)+3),') ',char(season_names(k))],'fontsize',16)
                xlabel('Surface area (m^2)','fontsize',16);
                if k == 2
                    ylabel('Normalized iceberg count','fontsize',16);
                end
            end
            %add the raw data
            yyaxis right; axr = gca; axr.XScale = 'linear';
            pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);
            maxn = max(log_n(:,:,l),[],2); minn = min(log_n(:,:,l),[],2);
            fill(gca,[log_v(~isnan(maxn)); flipud(log_v(~isnan(minn)))],[maxn(~isnan(maxn)); flipud(minn(~isnan(minn)))],...
                seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
            if contains(suffix,'zoom')
                set(axr,'ylim',[-6 3],'xlim',[1.5 5],'yticklabel',[]);
            else
                set(axr,'ylim',[-8 3],'xlim',[1.5 6.5],'yticklabel',[]);
            end
            clear maxn minn;
            axl.YAxis(2).Color = 'k';

            clear n1_mod v dv n n1 alpha c1 c2 c3 c4;
        end
        drawnow;

        %add the legend
        if l == 1 && k == 2
            seas_leg = legend(pf,season_names(2:4));
            set(seas_leg,'location','southoutside','orientation','horizontal')
            set(gca,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);
            set(axl,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);
            set(axr,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);
        end
        clear axl axr;
    end
    disp('moving on')
end
set(seas_leg,'position',[0.43 0.94 0.2,0.03]);
saveas(subdist_fig,[root_dir,'GrIS-melange_iceberg-distribution-variability_',suffix,'.png'],'png'); %save the plot
exportgraphics(subdist_fig,[root_dir,'GrIS-melange_iceberg-distribution-variability_',suffix,'.tif'],Resolution=600);

