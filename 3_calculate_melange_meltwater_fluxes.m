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
root_dir = '/Volumes/Jokulhaup_5T/Greenland-melange/'; %include trailing / in file name
% root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-melt/';
cd(root_dir);
EPSG_file = '/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj'; %generic EPSG3413 projection file

%create a month naming matrix to convert 3-letter months to numeric months
month_array = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};
years = 2011:1:2020; start_yr = years(1); end_yr = years(end);

%2-letter region flagging (based on alphabetical order of site folders)
site_abbrevs = ['ASG'];
site_names = [{'Alison'}];
region_flag = ['NW'];
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
disp('Converting iceberg size distribution textfiles for each date to date-specific csvs & a single csv');
convert_sizedistribution_txt_to_csv(root_dir);

disp('Melange size distributions for each site saved in *-melange-distributions.csv files');

%% Section 2: Export generic melange masks as shapefiles
disp('Export uncropped melange masks as shapefiles');
batch_export_melange_boundingbox_shapefiles(root_dir,EPSG_file);

%% Section 3: Estimate melange extent  
%Use velocity time series at points along the fjord (GRAB USING ITS_LIVE), 
%time-stamped melange masks, & the uncropped melange shapefile exported in
%Section 2


%% Section 4: Flag melange maps as partial or full for surface area estimates
close all;
disp('Checking DEMs to flag if they cover all or most of the melange');

for p = 1:size(site_abbrevs,1);
    disp(site_abbrevs(p,:));
    cd([root_dir,site_abbrevs(p,:)]);
    
    %create a dummy matrix to hold the melange extent flag
    melext(p).flag = []; %fill in with 2 when full melange in DEM, 1 when most, 0 otherwise
    
    %navigate to site folder and load data
    cd([root_dir,site_abbrevs(p,:)]);
    load([site_abbrevs(p,:),'-melange-masks.mat']); %load all the melange DEM masks
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
    
    %load each DEM & flag it
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
                melext(p).flag(1,k) = 2;
            case 'Most'
                melext(p).flag(1,k) = 1;
            case 'No'
                melext(p).flag(1,k) = 0;;
        end
        
        clear M extent melange; close all;
    end
    clear melange_mats melangemat_dates;
end
save([root_dir,'melange-extent-flags.mat'],'melext','-v7.3');

%% Section 3: Load iceberg size distributions & melange masks
close all;
disp('Bringing datasets together to estimate melange melt fluxes');
load([root_dir,'melange-extent-flags.mat']); %load the melange extent flags as needed
%  TESTING FOR A SINGLE SITE WITH SPARSE DATA... NEED TO CHECK AUTOMATION

%navigate to site folder and load data
for p = 1:size(site_abbrevs,1);
    disp(site_abbrevs(p,:));
    cd([root_dir,site_abbrevs(p,:)]);
%     load([site_abbrevs(p,:),'-melange-masks.mat']); %load all the melange DEM masks
    dists = readtable([site_abbrevs(p,:),'-melange-distributions.csv'],"VariableNamingRule","preserve"); %load all iceberg size distributions
    
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
        pl(1) = loglog(berg_A(~isnan(berg_fracs(:,j))),berg_fracs(~isnan(berg_fracs(:,j)),j),'-','color','k'); hold on;
    end
    pl(2) = loglog(berg_A,DJF_fracs(:,1),'-b','linewidth',2); hold on;
    pl(3) = loglog(berg_A,MAM_fracs(:,1),'-g','linewidth',2); hold on;
    pl(4) = loglog(berg_A,JJA_fracs(:,1),'-m','linewidth',2); hold on;
    pl(5) = loglog(berg_A,SON_fracs(:,1),'-r','linewidth',2); hold on;
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
    
    %estimate "typical" surface & submarged areas for each season
    DJF_SA = nanmedian(nansum(berg_nos(:,DJF(find(melext(p).flag(DJF)>=1))).*berg_A,1))
    MAM_SA = nanmedian(nansum(berg_nos(:,MAM(find(melext(p).flag(MAM)>=1))).*berg_A,1))
    JJA_SA = nanmedian(nansum(berg_nos(:,JJA(find(melext(p).flag(JJA)>=1))).*berg_A,1))
    SON_SA = nanmedian(nansum(berg_nos(:,SON(find(melext(p).flag(SON)>=1))).*berg_A,1))
%     SA = nanmedian(nansum(berg_nos.*berg_A,1)); %m^2
    DJF_subA = ((2*(900/1026)+1).*berg_A).*(DJF_fracs.*DJF_SA);
    MAM_subA = ((2*(900/1026)+1).*berg_A).*(MAM_fracs.*MAM_SA);
    JJA_subA = ((2*(900/1026)+1).*berg_A).*(JJA_fracs.*JJA_SA);
    SON_subA = ((2*(900/1026)+1).*berg_A).*(SON_fracs.*SON_SA);
    
    %multiply the surface area of each draft bin by the estimated best-fit 
    %melt rate for that draft
    berg_meltrate = interp1(draft_vector,meltcurve,berg_draft); berg_meltrate(end) = berg_meltrate(end-1);
    %binned meltwater flux (m^3/d)
    DJF_meltflux = berg_meltrate.*DJF_subA; 
    MAM_meltflux = berg_meltrate.*MAM_subA;
    JJA_meltflux = berg_meltrate.*JJA_subA;
    SON_meltflux = berg_meltrate.*SON_subA;
    
    %report seasonal meltwater fluxes in m^3/s
    fprintf('DJF typical meltwater flux = %4.0f m^3/s \n',nansum(DJF_meltflux)./86400);
    fprintf('MAM typical meltwater flux = %4.0f m^3/s \n',nansum(MAM_meltflux)./86400);
    fprintf('JJA typical meltwater flux = %4.0f m^3/s \n',nansum(JJA_meltflux)./86400);
    fprintf('SON typical meltwater flux = %4.0f m^3/s \n',nansum(SON_meltflux)./86400);
    
    %export the data to a csv table
    T=table(['DJF';'MAM';'JJA';'SON'],[DJF_SA; MAM_SA; JJA_SA; SON_SA],...
        [nansum(DJF_subA); nansum(MAM_subA); nansum(JJA_subA); nansum(SON_subA)],...
        [nansum(DJF_meltflux)./86400; nansum(MAM_meltflux)./86400; nansum(JJA_meltflux)./86400; nansum(SON_meltflux)./86400]);
    column_names = ["Months", "SurfaceArea", "SubmergedArea","MeltwaterFlux"];
    column_units = ["MMM", "m^2", "m^2","m^3/s"];
    T.Properties.VariableNames = column_names; T.Properties.VariableUnits = column_units;
    writetable(T,[root_dir,site_abbrevs(p,:),'/',site_abbrevs(p,:),'-area-meltflux-season-summary.csv']);
    
    clear dist* berg_* unique_dates refs inds melt* draft* DJF* MAM* JJA* SON* subA dVdt TF pl T;
    disp('... moving on to the next site');
end

