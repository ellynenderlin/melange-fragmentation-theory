%%% 3_calculate_melange_meltwater_fluxes: Use this code to combine
%%% estimated iceberg size distributions and melt rates to calculate fjord
%%% iceberg meltwater fluxes.

%% Section 0: Initialize (run every time)
clearvars; close all;
%add path for repo
addpath('/Users/ellynenderlin/Research/NSF_Greenland-Calving/iceberg-calving/DEMsizes_matlab-python');
%add paths for any supporting codes
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/LSQ_LUT_PIECEWISE');

%specify directories for required files ([root_dir,'/',site_abbrevs(i)])
basepath='/Volumes/Jokulhaup_5T/Greenland-melange/'; %this should be the overarching directory, with site-specific sub-directories
root_dir = basepath; output_dir = basepath;
cd(root_dir);
EPSG_file = '/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj'; %generic EPSG3413 projection file (available in group GitHub general-code)

%create a month naming matrix to convert 3-letter months to numeric months
month_array = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};
years = 2011:1:2020; start_yr = years(1); end_yr = years(end);

%2-letter region flagging (based on alphabetical order of site folders)
site_abbrevs = ['KBG'];
site_names = [{'Koge Bugt'}];
region_flag = ['SE'];
site_geog_order = [1]; %geographic order, counterclockwise from NW
% site_abbrevs = ['AG';'HH';'HM';'IB';'IG';'JI';'KB';'KL';'KO';'MD';'MG';'RI';'UM';'UN';'US';'ZI']; %alphabetical site directory list
% site_names = [{'Alison'},{'Helheim'},{'Ullip Sermia'},{'Salliarutsip Sermia'},{'Illulip Sermia'},...
%     {'Sermeq Kujalleq'},{'Koge Bugt'},{'Kangerlussuaq'},{'Kong Oscar'},{'Magga Dan'},...
%     {'Nigertiip Apusiia'},{'Kangilliup Sermia'},{'Umiammakku Sermia'},{'Upernavik North'},{'Upernavik South'},{'Zachariae Isstrom'}];
% region_flag = ['NW';'SE';'NW';'CW';'NW';'CW';'SE';'SE';'NW';'CE';'SE';'CW';'CW';'NW';'NW';'NE'];
% site_geog_order = [3;12;1;7;4;10;11;14;2;15;13;9;8;5;6;16]; %geographic order, counterclockwise from NW


%NOTE: data are added to F structure according to index specified by
%site_geog_order and plot locations and letters below are for those sorted data
%create a vector specifying subplot locations according to geography
site_plotlocs = ones(size(site_geog_order));
site_plotlocs(1:10) = 2*[1:1:10]'-1; %based on 10 sites along W coast
site_plotlocs(11:end-1) = [11:1:length(site_plotlocs)-1]' - 3*([11:1:length(site_plotlocs)-1]'-14);
site_plotlocs(end) = 2;
%letters for labels
site_plotletts = [{'a)'},{'b)'},{'c)'},{'d)'},{'e)'},{'f)'},{'g)'},{'h)'},{'i)'},{'j)'},{'k)'},{'l)'},{'m)'},{'n)'},{'o)'},{'p)'}];

%specify a regional colorramp
regions = unique(region_flag,'rows');
region_cmap = [171,217,233; 253,174,97; 44,123,182; 255,255,191; 215,25,28]/255; %same colorblind-friendly color scheme as plot_Greenland_iceberg_calving_fragmentation_fits.m

%specify annual colorramp
annual_cmap = cmocean('thermal',length([start_yr:1:end_yr]));
%specify seasonal colorramp
season_cmap = cmocean('curl',6); season_cmap = season_cmap(2:end-1,:);
%specify daily colorramp
% day_cmap = cmocean('phase',365);
accum_cmap = cmocean('ice',233); ablat_cmap = flipud(cmocean('solar',162));
day_cmap = [accum_cmap(103:end-10,:); ablat_cmap(11:end,:); accum_cmap(11:102,:)];


disp('Initialized iceberg meltwater flux code');
%% Section 1: Convert date-specific size distribution textfiles to a concatenated csv
disp('Converting iceberg size distribution textfiles for each date to date-specific csvs & a single csv...');
convert_sizedistribution_txt_to_csv(root_dir);

disp('Melange size distributions for each site saved in *iceberg-distributions.csv files');

%% Section 2: Export generic melange masks as shapefiles
disp('Exporting uncropped melange masks as shapefiles...');
batch_export_melange_boundingbox_shapefiles(root_dir,EPSG_file);
disp('DEM-based melange mask shapefiles created.');

disp('Create monthly delineations of the melange edge in GEEDiT');
disp(' & export as a shapefile to the "shapefiles" directory before running the next subsection');

%% Section 3: Manually estimate melange extent
close all;
disp('Create/compile manual melange extent estimates...');

%loop through each site
for p = 1:size(site_abbrevs,1);
    disp(site_abbrevs(p,:));
    
    %navigate to site folder and load melange mask data
    cd([root_dir,site_abbrevs(p,:)]);
    load([output_dir,site_abbrevs(p,:),'/',site_abbrevs(p,:),'-melange-masks.mat']); %load all the melange DEM masks
    for k = 1:length(melmask.dated); melA(k) = polyarea(melmask.dated(k).x,melmask.dated(k).y); end
    disp(['Typical fjord area (including open water): ',num2str(round(nanmean(melA))),'m^2']);
    
    %identify where the date is in the file name
    %account for different date location in file name depending on length of site abbreviation
    if length(site_abbrevs(p,:)) == 3
        matfile_daterefs = [5:12];
    elseif length(site_abbrevs(p,:)) == 2
        matfile_daterefs = [4:11];
    else
        error('Using a non-standard naming format! Switch to a 2- or 3-letter site abbreviation.');
    end
    
    %find DEMs
    melange_mats = dir([root_dir,site_abbrevs(p,:),'/DEMs/*_melange-DEMfilled.mat']);
    for i = 1:length(melange_mats); melangemat_dates(i,:) = melange_mats(i).name(matfile_daterefs); end
    disp([num2str(size(melange_mats,1)),' DEMs from ',melangemat_dates(1,1:4),'-',melangemat_dates(end,1:4)]);
    
    %load each DEM & flag it (fill in with 2 when full melange in DEM, 1 when most, 0 otherwise)
    for k = 1:length(melange_mats)
        disp(num2str(k));
        load([root_dir,site_abbrevs(p,:),'/DEMs/',melange_mats(k).name]);
        
        %extract melange elevations
        melange = M.DEM.z_filled;
        melange(isnan(M.DEM.z_filled)) = 0;
        melange(melange<0) = 0;
        melange(M.mask.DEM==0) = NaN; melange(M.mask.blunders==1) = NaN;
        disp([num2str(sum(~isnan(melange(melange>0)))*4),' m^2 with >0m elevations']);
        
        %plot the melange DEM
        close(gcf);
        figure; h = histogram(melange(~isnan(melange)),min(melange(~isnan(melange))):1:max(melange(~isnan(melange)))); %use the elevation histogram to set colormap range
        cmax = h.BinEdges(find(cumsum(h.Values)>=0.98*sum(h.Values),1,'first')+1);
        clear h; close(gcf);
        figDEM = figure; set(figDEM,'position',[550 100 1000 500]);
        imagesc(M.DEM.x,M.DEM.y,melange); axis xy equal;
        melange_cmap = cmocean('thermal',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap);
        hold on; DEMax = gca;
        set(gca,'clim',[0 16]);  %set(gca,'clim',[0 cmax]);
        cbar = colorbar; cbar.Label.String  = 'elevation (m a.s.l.)';
        set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]); xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        title([melangemat_dates(k,1:4),'/',melangemat_dates(k,5:6),'/',melangemat_dates(k,7:8)],'fontsize',16); grid on; drawnow;
        
        %fill in the extent matrix
        extent = questdlg('Does the melange extent look complete?',...
            'Extent Check','Yes','Most','No','No');
        switch extent
            case 'Yes'
                melmask.dated(k).DEMext_flag = 2;
            case 'Most'
                melmask.dated(k).DEMext_flag = 1;
            case 'No'
                melmask.dated(k).DEMext_flag = 0;;
        end
        
        clear M extent melange; close all;
    end
    clear melange_mats melangemat_dates melange;
    save([output_dir,site_abbrevs(p,:),'/',site_abbrevs(p,:),'-melange-masks.mat'],'melmask','-v7.3');
    
    
    %bring in manually-delineated melange extents from GEEDiT for comparison
    %(if they exist)
    cd([root_dir,site_abbrevs(p,:),'/shapefiles/']);
    sfile = dir('*melange-delineations*.shp'); %find the melange shapefile & load it (must have "melange-delineations" in name)
    if ~isempty(sfile) %if the directory containing GEEDiT shapefiles exists, grab the data
        disp('Incorporating manual image-based delineations into melange extent estimates');
        S = shaperead(sfile(1).name);
        for k = 1:size(S,1)
            if length(S(k).Date) == 10
                melext(k).date = [S(k).Date(1:4),S(k).Date(6:7),S(k).Date(9:10)];
                ext_decidate(k,:) = convert_to_decimaldate(melext(k).date);
            else
                error('Unfamiliar date formatting from GEEDiT');
            end
            melext(k).lat = S(k).Y; melext(k).lon = S(k).X;
            [melext(k).x,melext(k).y] = wgs2ps(melext(k).lon,melext(k).lat);
        end
        clear S;
        save([output_dir,site_abbrevs(p,:),'/',site_abbrevs(p,:),'-melange-delineations.mat'],'melext','-v7.3');
        
        %crop the time-stamped melange outlines using the delineations
%         figure;
        for k = 1:length(melmask.dated);
            %idenfity the delineation with the closest time stamp
            DEM_decidate(k,:) = convert_to_decimaldate(melmask.dated(k).datestring);
            date_diff = DEM_decidate(k,:) - ext_decidate;
            if min(abs(date_diff)) > 1; %if there is more than a year between the most recent terminus delineation & DEM
                clear date_diff; date_diff = (DEM_decidate(k,:)-floor(DEM_decidate(k,:))) - (ext_decidate-floor(ext_decidate));
                [sorted_dates,inds] = sort(abs(date_diff));
                date_ref = inds(1:2); %select the two dates that are seasonally closest
                clear sorted_dates inds;
            else
                date_ref = find(abs(date_diff)==min(abs(date_diff)),1,'first');
            end
            
            %calculate the area of each polygon that approximates the 
            for j = 1:length(date_ref)
                out_intercept = []; out_interceptx = []; out_intercepty = [];
                for i = 1:length(melmask.dated(k).x)-1
                    [xi,yi] = polyxpoly(melmask.dated(k).x(i:i+1),melmask.dated(k).y(i:i+1),melext(date_ref(j)).x,melext(date_ref(j)).y); %find the intersections of the terminus trace with the melange outline
                    if ~isempty(xi)
                        out_intercept = [out_intercept i]; out_interceptx = [out_interceptx xi]; out_intercepty = [out_intercepty yi];
                    end
                    clear xi yi;
                end
                
                %create a new melange polygon cropped on both ends
                x_mel = [out_interceptx(find(out_intercept==min(out_intercept))); melmask.dated(k).x(min(out_intercept)+1:max(out_intercept)-1); out_interceptx(find(out_intercept==max(out_intercept)));out_interceptx(find(out_intercept==min(out_intercept)))];
                y_mel = [out_intercepty(find(out_intercept==min(out_intercept))); melmask.dated(k).y(min(out_intercept)+1:max(out_intercept)-1); out_intercepty(find(out_intercept==max(out_intercept)));out_intercepty(find(out_intercept==min(out_intercept)))];
                melpoly = polyshape(x_mel,y_mel); SA(j) = area(melpoly);
                clear x_mel y_mel melpoly;
            end
            clear date_diff date_ref;
            melmask.dated(k).IMext_area = median(SA); melmask.dated(k).IMext_areaMAD = mad(SA,1); clear SA;
%             plot(DEM_decidate(k,:),melmask.dated(k).IMext_area,'xk'); hold on; drawnow;
            clear *intercept*;
        end
        save([output_dir,site_abbrevs(p,:),'/',site_abbrevs(p,:),'-melange-masks.mat'],'melmask','-v7.3');
    else
        disp('Delineate seaward melange edge every month from 2016-present in GEEDiT!');
    end
    clear sfile;
    disp('Surface areas from manual delineations added to melmask structure');
    
end
clear melmask melext;
disp('Move on to the last section!');

%% Section 4: Load iceberg size distributions & melange extent estimates
close all;
disp('Bringing datasets together to estimate melange melt fluxes');
%  TESTING FOR A SINGLE SITE... NEED TO CHECK AUTOMATION

%seasonal colormap
season_cmap = [233,163,201; 161,215,106; 77,146,33; 197,27,125]./255; %green-pink from colorbrewer

%specify which dataset you want to use to estimate melange surface areas
Asource = questdlg('Specify if you want to use the velocity coherence, manual melange extents, or both.',...
    'Melange Extent Data','Manual','Velocity','Both','Manual');

%navigate to site folder and load data
for p = 1:size(site_abbrevs,1);
    disp(site_abbrevs(p,:));
    cd([root_dir,site_abbrevs(p,:)]);
    dists = readtable([site_abbrevs(p,:),'-iceberg-distributions.csv'],"VariableNamingRule","preserve"); %load all iceberg size distributions
    load([output_dir,site_abbrevs(p,:),'/',site_abbrevs(p,:),'-melange-masks.mat']); %load all dated melange datasets
    load([output_dir,site_abbrevs(p,:),'/',site_abbrevs(p,:),'-melange-delineations.mat']); %load manual melange delineations
    
    %use iceberg size distributions to estimate the submerged area of the
    %melange in each season (DJF, MAM, JJA, SON)
    dist_hdrs = dists.Properties.VariableNames;
    berg_A = dists.SurfaceArea_mean; berg_dA = dists.SurfaceArea_range; %circular surface areas (m^2)
    for j = 1:length(dist_hdrs)-2
        berg_datestr(j,:) = char(dist_hdrs(j+2));
        berg_yr(j,:) = str2num(berg_datestr(j,1:4));
        berg_mo(j,:) = str2num(berg_datestr(j,6:7));
        berg_decidate(j,:) = convert_to_decimaldate([berg_datestr(j,1:4),berg_datestr(j,6:7),berg_datestr(j,9:10)]);
    end
    berg_draft = (900/1026)*sqrt(berg_A/pi());
    %create a table of iceberg distribution numbers & normalize by area
    berg_nos = table2array(dists(:,3:end));
    berg_fracs = berg_nos./(nansum(berg_nos.*berg_A,1));
    clear dists;
    
    %identify the seasons for the observations
    DJF = find(berg_mo==12 | berg_mo == 1 | berg_mo == 2);
    MAM = find(berg_mo==3 | berg_mo == 4 | berg_mo == 5);
    JJA = find(berg_mo==6 | berg_mo == 7 | berg_mo == 8);
    SON = find(berg_mo==9 | berg_mo == 10 | berg_mo == 11);
    
    %solve for median fractional distributions & plot
    DJF_fracs = nanmedian(berg_fracs(:,DJF),2);
    MAM_fracs = nanmedian(berg_fracs(:,MAM),2);
    JJA_fracs = nanmedian(berg_fracs(:,JJA),2);
    SON_fracs = nanmedian(berg_fracs(:,SON),2);
    figure; set(gcf,'position',[50 50 600 600]);
    for j = 1:size(berg_fracs,2)
        pl(1) = loglog(berg_A(~isnan(berg_fracs(:,j))),berg_fracs(~isnan(berg_fracs(:,j)),j),'-','color',[0.5 0.5 0.5]); hold on;
    end
    pl(2) = loglog(berg_A,DJF_fracs(:,1),'-','linewidth',2,'color',season_cmap(1,:)); hold on;
    pl(3) = loglog(berg_A,MAM_fracs(:,1),'-','linewidth',2,'color',season_cmap(2,:)); hold on;
    pl(4) = loglog(berg_A,JJA_fracs(:,1),'-','linewidth',2,'color',season_cmap(3,:)); hold on;
    pl(5) = loglog(berg_A,SON_fracs(:,1),'-','linewidth',2,'color',season_cmap(4,:)); hold on;
    grid on; set(gca,'fontsize',16);
    leg=legend(pl,'all','DJF','MAM','JJA','SON');
    xlabel('Iceberg surface area (m^2)','fontsize',16);
    ylabel('Count normalized by total area','fontsize',16);
    
    %load all meltrates & concatenate
    meltdates = []; draft = []; subA = []; dVdt = [];
    meltfiles = dir([root_dir,site_abbrevs(p,:),'/meltrates/*.csv']);
    for j = 1:length(meltfiles);
        melts = readtable([meltfiles(j).folder,'/',meltfiles(j).name]);
        draft = [draft; melts.MedianDraft_mean];  %m b.s.l.
        subA = [subA; melts.SubmergedArea_mean]; %m^2
        dVdt = [dVdt; melts.VolumeChangeRate]; %m^3/d
        meltdates = [meltdates; [repmat(string(meltfiles(j).name(length(site_abbrevs(p,:))+2:length(site_abbrevs(p,:))+9)),size(melts.MedianDraft_mean)),...
            repmat(string(meltfiles(j).name(length(site_abbrevs(p,:))+11:length(site_abbrevs(p,:))+18)),size(melts.MedianDraft_mean))]];
        clear melts;
    end
    meltrate = dVdt./subA; %m/d
    
    %visually check that the meltrate data look reasonable
    [unique_dates,refs,inds] = unique(meltdates,'rows');
    figure; set(gcf,'position',[650 50 600 600]);
    cmap = cmocean('matter',size(refs,1)); 
    for k = 1:length(inds); colors(k,:) = cmap(inds(k),:); end
    scatter(draft,meltrate,24,colors,'filled'); hold on; grid on;
    
    %fit a complex curve to the meltrate vs draft data to parameterize melt
    draft_vector = [0:1:ceil(max(berg_draft))];
    ft = fittype('poly4'); %try a 4th order polynomial to capture expected melt variability with draft
    opts = fitoptions('Method','LinearLeastSquares' );
    opts.Lower = [-Inf -Inf -Inf -Inf 0]; opts.Upper = [Inf Inf Inf Inf 0]; %force intercept through 0,0
    [fitresult,gof] = fit(draft,meltrate,ft,opts); meltcurve = feval(fitresult,draft_vector);
    plot(draft_vector(1:round(max(draft))+1),meltcurve(1:round(max(draft))+1),'-k'); hold on; %add fit curve to scatterplot of data
    %find the local maximum for deep icebergs and fix that value for all deeper icebergs
    TF = islocalmax(meltcurve(1:round(max(draft))+1)); %only find maximum over the oberved melt draft range 
    TF(round(max(draft))+2:length(draft_vector)) = 0;
    draft_meltmax = draft_vector(find(TF==1,1,'last'));
    meltrate_meltmax = meltcurve(find(TF==1,1,'last'));
    meltcurve(find(TF==1,1,'last'):end) = meltcurve(find(TF==1,1,'last'));
    plot(draft_vector(1:round(max(draft))+1),meltcurve(1:round(max(draft))+1),'--k','linewidth',2); hold on; %add adjusted fit curve to scatterplot of data
    set(gca,'fontsize',16);
    xlabel('Iceberg draft (m b.s.l.)','fontsize',16);
    ylabel('Melt rate (m/d)','fontsize',16);
    
    %estimate "typical" surface & submerged areas for each season
    switch Asource
        case 'Manual'
            %area estimates
            disp('DEM and image-based area estimates');
            %DEM-based
            for k = 1:length(melmask.dated); DEMext_flag(k) = melmask.dated(k).DEMext_flag; end
            DJF_SA(1) = nanmedian(nansum(berg_nos(:,DJF(find(DEMext_flag(DJF)>=1))).*berg_A,1));
            DJF_SAmad(1) = mad(nansum(berg_nos(:,DJF(find(DEMext_flag(DJF)>=1))).*berg_A,1),1);
            MAM_SA(1) = nanmedian(nansum(berg_nos(:,MAM(find(DEMext_flag(MAM)>=1))).*berg_A,1));
            MAM_SAmad(1) = mad(nansum(berg_nos(:,MAM(find(DEMext_flag(MAM)>=1))).*berg_A,1),1);
            JJA_SA(1) = nanmedian(nansum(berg_nos(:,JJA(find(DEMext_flag(JJA)>=1))).*berg_A,1));
            JJA_SAmad(1) = mad(nansum(berg_nos(:,JJA(find(DEMext_flag(JJA)>=1))).*berg_A,1),1);
            SON_SA(1) = nanmedian(nansum(berg_nos(:,SON(find(DEMext_flag(SON)>=1))).*berg_A,1));
            SON_SAmad(1) = mad(nansum(berg_nos(:,SON(find(DEMext_flag(SON)>=1))).*berg_A,1),1);
            %image-delineated extents if available
            if isfield(melmask.dated(k),'IMext_area')
                for k = 1:length(melmask.dated); SA(k) = melmask.dated(k).IMext_area; SAmad(k) = melmask.dated(k).IMext_areaMAD; end
                DJF_SA(2) = nanmedian(SA(DJF)); DJF_SAmad(2) = max([SAmad,mad(SA(DJF),1)]); 
                MAM_SA(2) = nanmedian(SA(MAM)); MAM_SAmad(2) = max([SAmad,mad(SA(MAM),1)]); 
                JJA_SA(2) = nanmedian(SA(JJA)); JJA_SAmad(2) = max([SAmad,mad(SA(JJA),1)]); 
                SON_SA(2) = nanmedian(SA(SON)); SON_SAmad(2) = max([SAmad,mad(SA(SON),1)]); 
                clear SA SAmad:
            end
            clear DEMext_flag;
            %velocity-based
            DJF_SA(3) = NaN; MAM_SA(3) = NaN; JJA_SA(3) = NaN; SON_SA(3) = NaN; 
            
            %calculate the average submerged area from the DEMs & manual
            %delineations (DEMs may include sea ice, manual delineations
            %are subjective when the ice coverage is sparse)
            DJF_subA = ((2*(900/1026)+1).*berg_A).*(DJF_fracs.*nanmean(DJF_SA(1:2)));
            DJF_subA_range = ((2*(900/1026)+1).*berg_A).*(DJF_fracs.*[(nanmean(DJF_SA(1:2))-max(DJF_SAmad(1:2))), (nanmean(DJF_SA(1:2))+max(DJF_SAmad(1:2)))]);
            MAM_subA = ((2*(900/1026)+1).*berg_A).*(MAM_fracs.*nanmean(MAM_SA(1:2)));
            MAM_subA_range = ((2*(900/1026)+1).*berg_A).*(MAM_fracs.*[(nanmean(MAM_SA(1:2))-max(MAM_SAmad(1:2))), (nanmean(MAM_SA(1:2))+max(MAM_SAmad(1:2)))]);
            JJA_subA = ((2*(900/1026)+1).*berg_A).*(JJA_fracs.*nanmean(JJA_SA(1:2)));
            JJA_subA_range = ((2*(900/1026)+1).*berg_A).*(JJA_fracs.*[(nanmean(JJA_SA(1:2))-max(JJA_SAmad(1:2))), (nanmean(JJA_SA(1:2))+max(JJA_SAmad(1:2)))]);
            SON_subA = ((2*(900/1026)+1).*berg_A).*(SON_fracs.*nanmean(SON_SA(1:2)));
            SON_subA_range = ((2*(900/1026)+1).*berg_A).*(SON_fracs.*[(nanmean(SON_SA(1:2))-max(SON_SAmad(1:2))), (nanmean(SON_SA(1:2))+max(SON_SAmad(1:2)))]);
            
            %report seasonal surface areas in km^2
            fprintf('DJF median (MAD) extent = %3.0f (%3.0f) km^2 \n',nanmean(DJF_SA(1:2))./10^6,max(DJF_SAmad(1:2))./10^6);
            fprintf('MAM median (MAD) extent = %3.0f (%3.0f) km^2 \n',nanmean(MAM_SA(1:2))./10^6,max(MAM_SAmad(1:2))./10^6);
            fprintf('JJA median (MAD) extent = %3.0f (%3.0f) km^2 \n',nanmean(JJA_SA(1:2))./10^6,max(JJA_SAmad(1:2))./10^6);
            fprintf('SON median (MAD) extent = %3.0f (%3.0f) km^2 \n',nanmean(SON_SA(1:2))./10^6,max(SON_SAmad(1:2))./10^6);
            
        case 'Velocity'
            %Use velocity time series at points along the fjord (from ITS_LIVE),
            %time-stamped melange masks, & the uncropped melange shapefile exported in Section 2
            disp('Estimating melange extent from ITS_LIVE velocity time series...');
            %if the velocities directory exists, go into it and check for files
            if exist('velocities') == 7
                cd('velocities/');
                vels = dir('*.csv');
                
                %specify the path and filename for uncropped fjord outline & dated melange outlines
                outline_dir = [output_dir,site_abbrevs(p,:),'/']; outline_file = [site_abbrevs(p,:),'-melange-masks.mat'];
                
                %loop through the velocity data and add melange areas based on
                %dates with velocity (indicating coherent melange) to the melmask
                %structure
                if ~isempty(vels)
                    [melmask,velSA,velSAmad] = analyze_melange_velocity_coherence(site_abbrevs(p,:),root_dir,outline_dir,outline_file,melmask);
                    DJF_SA(3) = velSA(1); MAM_SA(3) = velSA(2); JJA_SA(3) = velSA(3); SON_SA(3) = velSA(4);
                    DJF_SAmad(3) = velSAmad(1); MAM_SAmad(3) = velSAmad(2); JJA_SAmad(3) = velSAmad(3); SON_SAmad(3) = velSAmad(4);
                else
%                     DJF_SA(3) = NaN; MAM_SA(3) = NaN; JJA_SA(3) = NaN; SON_SA(3) = NaN;
%                     DJF_SAmad(3) = NaN; MAM_SAmad(3) = NaN; JJA_SAmad(3) = NaN; SON_SAmad(3) = NaN;
                    error('No velocity csv files... grab time series from points using ITS_LIVE portal & add to velocities directory');
                end
                clear outline_* vels;
            else
                error('No velocities directory and csv files... grab time series from points using ITS_LIVE portal');
            end

            %calculate the average submerged area from the velocity coherence extents
            DJF_subA = ((2*(900/1026)+1).*berg_A).*(DJF_fracs.*DJF_SA(3));
            DJF_subA_range = ((2*(900/1026)+1).*berg_A).*(DJF_fracs.*[(DJF_SA(3)-DJF_SAmad(3)), (DJF_SA(3)+DJF_SAmad(3))]);
            MAM_subA = ((2*(900/1026)+1).*berg_A).*(MAM_fracs.*MAM_SA(3));
            MAM_subA_range = ((2*(900/1026)+1).*berg_A).*(MAM_fracs.*[(MAM_SA(3)-MAM_SAmad(3)), (MAM_SA(3)+MAM_SAmad(3))]);
            JJA_subA = ((2*(900/1026)+1).*berg_A).*(JJA_fracs.*JJA_SA(3));
            JJA_subA_range = ((2*(900/1026)+1).*berg_A).*(JJA_fracs.*[(JJA_SA(3)-JJA_SAmad(3)), (JJA_SA(3)+JJA_SAmad(3))]);
            SON_subA = ((2*(900/1026)+1).*berg_A).*(SON_fracs.*SON_SA(3));
            SON_subA_range = ((2*(900/1026)+1).*berg_A).*(SON_fracs.*[(SON_SA(3)-SON_SAmad(3)), (SON_SA(3)+SON_SAmad(3))]);
            
            %report seasonal surface areas in km^2
            fprintf('DJF median (MAD) extent = %3.0f (%3.0f) km^2 \n',DJF_SA(3)./10^6,DJF_SAmad(3)./10^6);
            fprintf('MAM median (MAD) extent = %3.0f (%3.0f) km^2 \n',MAM_SA(3)./10^6,MAM_SAmad(3)./10^6);
            fprintf('JJA median (MAD) extent = %3.0f (%3.0f) km^2 \n',JJA_SA(3)./10^6,JJA_SAmad(3)./10^6);
            fprintf('SON median (MAD) extent = %3.0f (%3.0f) km^2 \n',SON_SA(3)./10^6,SON_SAmad(3)./10^6);
            
        case 'Both'
            %area estimates
            disp('Velocity, DEM, & image-based area estimates');
            %DEM-based
            for k = 1:length(melmask.dated); DEMext_flag(k) = melmask.dated(k).DEMext_flag; end
            DJF_SA(1) = nanmedian(nansum(berg_nos(:,DJF(find(DEMext_flag(DJF)>=1))).*berg_A,1));
            DJF_SAmad(1) = mad(nansum(berg_nos(:,DJF(find(DEMext_flag(DJF)>=1))).*berg_A,1),1);
            MAM_SA(1) = nanmedian(nansum(berg_nos(:,MAM(find(DEMext_flag(MAM)>=1))).*berg_A,1));
            MAM_SAmad(1) = mad(nansum(berg_nos(:,MAM(find(DEMext_flag(MAM)>=1))).*berg_A,1),1);
            JJA_SA(1) = nanmedian(nansum(berg_nos(:,JJA(find(DEMext_flag(JJA)>=1))).*berg_A,1));
            JJA_SAmad(1) = mad(nansum(berg_nos(:,JJA(find(DEMext_flag(JJA)>=1))).*berg_A,1),1);
            SON_SA(1) = nanmedian(nansum(berg_nos(:,SON(find(DEMext_flag(SON)>=1))).*berg_A,1));
            SON_SAmad(1) = mad(nansum(berg_nos(:,SON(find(DEMext_flag(SON)>=1))).*berg_A,1),1);
            %image-delineated extents if available
            if isfield(melmask.dated(k),'IMext_area')
                for k = 1:length(melmask.dated); SA(k) = melmask.dated(k).IMext_area; SAmad(k) = melmask.dated(k).IMext_areaMAD; end
                DJF_SA(2) = nanmedian(SA(DJF)); DJF_SAmad(2) = max([SAmad,mad(SA(DJF),1)]); 
                MAM_SA(2) = nanmedian(SA(MAM)); MAM_SAmad(2) = max([SAmad,mad(SA(MAM),1)]); 
                JJA_SA(2) = nanmedian(SA(JJA)); JJA_SAmad(2) = max([SAmad,mad(SA(JJA),1)]); 
                SON_SA(2) = nanmedian(SA(SON)); SON_SAmad(2) = max([SAmad,mad(SA(SON),1)]); 
                clear SA SAmad:
            end
            clear DEMext_flag;
            %velocity-based
            disp('Estimating melange extent from ITS_LIVE velocity time series...');
            %if the velocities directory exists, go into it and check for files
            if exist('velocities') == 7
                cd('velocities/');
                vels = dir('*.csv');
                
                %specify the path and filename for uncropped fjord outline & dated melange outlines
                outline_dir = [output_dir,site_abbrevs(p,:),'/']; outline_file = [site_abbrevs(p,:),'-melange-masks.mat'];
                
                %loop through the velocity data and add melange areas based on
                %dates with velocity (indicating coherent melange) to the melmask
                %structure
                if ~isempty(vels)
                    [melmask,velSA,velSAmad] = analyze_melange_velocity_coherence(site_abbrevs(p,:),root_dir,outline_dir,outline_file,melmask);
                    DJF_SA(3) = velSA(1); MAM_SA(3) = velSA(2); JJA_SA(3) = velSA(3); SON_SA(3) = velSA(4);
                    DJF_SAmad(3) = velSAmad(1); MAM_SAmad(3) = velSAmad(2); JJA_SAmad(3) = velSAmad(3); SON_SAmad(3) = velSAmad(4);
                else
%                     DJF_SA(3) = NaN; MAM_SA(3) = NaN; JJA_SA(3) = NaN; SON_SA(3) = NaN;
%                     DJF_SAmad(3) = NaN; MAM_SAmad(3) = NaN; JJA_SAmad(3) = NaN; SON_SAmad(3) = NaN;
                    error('No velocity csv files... grab time series from points using ITS_LIVE portal & add to velocities directory');
                end
                clear outline_* vels;
            else
                error('No velocities directory and csv files... grab time series from points using ITS_LIVE portal');
            end
            
            %calculate the average submerged area from the DEMs & manual
            %delineations (DEMs may include sea ice, manual delineations
            %are subjective when the ice coverage is sparse)
            DJF_subA = ((2*(900/1026)+1).*berg_A).*(DJF_fracs.*nanmean(DJF_SA(1:2)));
            DJF_subA_range = ((2*(900/1026)+1).*berg_A).*(DJF_fracs.*[(nanmean(DJF_SA(1:2))-max(DJF_SAmad(1:2))), (nanmean(DJF_SA(1:2))+max(DJF_SAmad(1:2)))]);
            MAM_subA = ((2*(900/1026)+1).*berg_A).*(MAM_fracs.*nanmean(MAM_SA(1:2)));
            MAM_subA_range = ((2*(900/1026)+1).*berg_A).*(MAM_fracs.*[(nanmean(MAM_SA(1:2))-max(MAM_SAmad(1:2))), (nanmean(MAM_SA(1:2))+max(MAM_SAmad(1:2)))]);
            JJA_subA = ((2*(900/1026)+1).*berg_A).*(JJA_fracs.*nanmean(JJA_SA(1:2)));
            JJA_subA_range = ((2*(900/1026)+1).*berg_A).*(JJA_fracs.*[(nanmean(JJA_SA(1:2))-max(JJA_SAmad(1:2))), (nanmean(JJA_SA(1:2))+max(JJA_SAmad(1:2)))]);
            SON_subA = ((2*(900/1026)+1).*berg_A).*(SON_fracs.*nanmean(SON_SA(1:2)));
            SON_subA_range = ((2*(900/1026)+1).*berg_A).*(SON_fracs.*[(nanmean(SON_SA(1:2))-max(SON_SAmad(1:2))), (nanmean(SON_SA(1:2))+max(SON_SAmad(1:2)))]);
            
            %report seasonal surface areas in km^2
            fprintf('DJF median (MAD) extent = %3.0f (%3.0f) km^2 \n',nanmean(DJF_SA(1:2))./10^6,max(DJF_SAmad(1:2))./10^6);
            fprintf('MAM median (MAD) extent = %3.0f (%3.0f) km^2 \n',nanmean(MAM_SA(1:2))./10^6,max(MAM_SAmad(1:2))./10^6);
            fprintf('JJA median (MAD) extent = %3.0f (%3.0f) km^2 \n',nanmean(JJA_SA(1:2))./10^6,max(JJA_SAmad(1:2))./10^6);
            fprintf('SON median (MAD) extent = %3.0f (%3.0f) km^2 \n',nanmean(SON_SA(1:2))./10^6,max(SON_SAmad(1:2))./10^6);
            
    end
    
    %multiply the surface area of each draft bin by the estimated best-fit 
    %melt rate for that draft
    berg_meltrate = interp1(draft_vector,meltcurve,berg_draft); berg_meltrate(end) = berg_meltrate(end-1);
    %binned meltwater flux (m^3/d)
    DJF_meltflux = berg_meltrate.*DJF_subA; DJF_meltflux_range = berg_meltrate.*DJF_subA_range; 
    MAM_meltflux = berg_meltrate.*MAM_subA; MAM_meltflux_range = berg_meltrate.*MAM_subA_range; 
    JJA_meltflux = berg_meltrate.*JJA_subA; JJA_meltflux_range = berg_meltrate.*JJA_subA_range; 
    SON_meltflux = berg_meltrate.*SON_subA; SON_meltflux_range = berg_meltrate.*SON_subA_range; 
    
    %report seasonal meltwater fluxes in m^3/s
    fprintf('DJF typical meltwater flux = %3.0f (%3.0f, %3.0f) m^3/s \n',nansum(DJF_meltflux)./86400,min(nansum(DJF_meltflux_range)./86400),max(nansum(DJF_meltflux_range)./86400));
    fprintf('MAM typical meltwater flux = %3.0f (%3.0f, %3.0f) m^3/s \n',nansum(MAM_meltflux)./86400,min(nansum(MAM_meltflux_range)./86400),max(nansum(MAM_meltflux_range)./86400));
    fprintf('JJA typical meltwater flux = %3.0f (%3.0f, %3.0f) m^3/s \n',nansum(JJA_meltflux)./86400,min(nansum(JJA_meltflux_range)./86400),max(nansum(JJA_meltflux_range)./86400));
    fprintf('SON typical meltwater flux = %3.0f (%3.0f, %3.0f) m^3/s \n',nansum(SON_meltflux)./86400,min(nansum(SON_meltflux_range)./86400),max(nansum(SON_meltflux_range)./86400));
    
    %export the data to a csv table
    T=table(['DJF';'MAM';'JJA';'SON'],[DJF_SA; MAM_SA; JJA_SA; SON_SA],...
        [nansum(DJF_subA); nansum(MAM_subA); nansum(JJA_subA); nansum(SON_subA)],...
        [nansum(DJF_meltflux)./86400; nansum(MAM_meltflux)./86400; nansum(JJA_meltflux)./86400; nansum(SON_meltflux)./86400]);
    column_names = ["Months", "SurfaceArea", "SubmergedArea","MeltwaterFlux"];
    column_units = ["MMM", "m^2", "m^2","m^3/s"];
    T.Properties.VariableNames = column_names; T.Properties.VariableUnits = column_units;
    writetable(T,[root_dir,site_abbrevs(p,:),'/',site_abbrevs(p,:),'-area-meltflux-season-summary.csv']);
    
    clear dist* berg_* unique_dates refs inds melt* draft* DJF* MAM* JJA* SON* *SA *subA dVdt ft gof opts TF pl T melmask melext;
    disp('... moving on to the next site');
end

