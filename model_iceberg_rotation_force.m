% calculate theoretical buttressing needed for icebergs trying to calve
% from termini & with varying thickness and rotation angles (keep aspect
% ratio fixed for simplicity)
clearvars; close all;

%specify variable ranges, paths, and constants
output_dir ='/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';
% H_range = [100:50:1000]; %range of iceberg thicknesses
ang_range = [1:1:5]; %range of angles from vertical
E = 0.42; %width-to-thickness ratio
rho_i = 900; rho_w = 1026; %densities
H_maxes = []; %dummy variables to hold maximum iceberg thickness estimates
H_range = []; %dummy variables to hold quartiles & max terminus thickness estimates

%read in the melange data
load([output_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']);

%read in BedMachine Greenland v5 (elevations referenced to the geoid)
cd('/Users/ellynenderlin/Research/miscellaneous/elevations/BedMachine/');
bx = ncread('BedMachineGreenland-v5.nc','x');
by = ncread('BedMachineGreenland-v5.nc','y');
bz = ncread('BedMachineGreenland-v5.nc','bed'); %bed elevations
bze = ncread('BedMachineGreenland-v5.nc','errbed'); %bed elevations
bzs = ncread('BedMachineGreenland-v5.nc','source'); %data interp (0-2 = direct & mass conservation)
bzi = ncread('BedMachineGreenland-v5.nc','dataid'); %data source (>2 = seismic or multibeam bathy)
bg = ncread('BedMachineGreenland-v5.nc','geoid'); %geoid elevations relative to WGS84 ellipsoid


% loop through each site and grab floatation-based thickness estimates
cd(output_dir);
contents = dir(output_dir);
folders = contents([contents.isdir]);
validFolders = folders(~ismember({folders.name}, {'.', '..'}));
folderNames = {validFolders.name};

for i = 1:length(folderNames)
    site_abbrev = char(folderNames(i));
    cd(site_abbrev);
    disp(site_abbrev);

    %account for different date location in file name depending on length of site abbreviation
    if length(site_abbrev) == 3
        matfile_daterefs = [5:12];
    elseif length(site_abbrev) == 2
        matfile_daterefs = [4:11];
    else
        error('Using a non-standard naming format! Switch to a 2- or 3-letter site abbreviation.');
    end

    %load supporting data to identify the terminus position along the centerline
    load([site_abbrev,'-melange-masks.mat']); %created using create_melange_masks
    mask_dates = '';
    for p = 1:length(melmask.dated)
        mask_dates(p,:) = melmask.dated(p).datestring;
    end
    cd('shapefiles/')
    S = shaperead([site_abbrev,'_centerline.shp']);

    %subset the bedmap to the melange area
    BB = [min(melmask.uncropped.x),min(melmask.uncropped.y); min(melmask.uncropped.x),max(melmask.uncropped.y);...
        max(melmask.uncropped.x),max(melmask.uncropped.y); max(melmask.uncropped.x),min(melmask.uncropped.y)];
    bed_xind = [find(bx<=min(BB(:,1)),1,'last'),find(bx>=max(BB(:,1)),1,'first')];
    bed_yind = [find(by>=max(BB(:,2)),1,'last'),find(by<=min(BB(:,2)),1,'first')];
    bed_xsub = bx(min(bed_xind):max(bed_xind)); bed_ysub = by(min(bed_yind):max(bed_yind));
    bed_zsub = bz(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))'; %rotate so rows are y, columns are x
    bed_zesub = bze(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))';
    bed_sourcesub = bzs(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))';
    bed_idsub = bzi(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))';
    bed_gsub = bg(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))'; %rotate so rows are y, columns are x
    
    %filter wonky data from bad interpolation & plot
    bed_sourcefilt = bed_zsub; bed_sourcefilt(bed_sourcesub>2) = NaN;
    bed_errfilt = bed_zsub; bed_errfilt(bed_zesub>=200) = NaN;
    % figure; imagesc(bed_xsub,bed_ysub,bed_zsub); axis xy equal; colormap jet; colorbar; hold on; %plot to check subsetting
    % for j = 1:length(melmask.dated)
    %     plot(melmask.dated(j).x,melmask.dated(j).y,'k'); hold on;
    % end
    % drawnow;

    %loop through the dates to grab terminus thickness estimates
    cd ('../DEMs/');
    DEM_mats = dir([site_abbrev,'*_melange-DEMfilled.mat']); DEM_dates = ''; %DEMs
    Hmax = []; Hiqr = [];
    for p = 1:length(DEM_mats)
        DEM_dates(p,:) = DEM_mats(p).name(matfile_daterefs);
        % disp(DEM_dates(p,:))
        load([output_dir,'/',site_abbrev,'/DEMs/',DEM_mats(p).name]);
        %     disp('Data loaded.');

        %fill data holes and masked regions with NaNs
        melange = M.DEM.z_filled;
        melange(melange<0) = NaN; %melange(melange>Hmax/(rho_i/(rho_sw-rho_i))) = NaN;
        melange(M.mask.DEM==0) = NaN;

        %interpolate to the same grid
        [z_xgrid,z_ygrid] = meshgrid(M.DEM.x,M.DEM.y);
        [b_xgrid,b_ygrid] = meshgrid(bed_xsub,bed_ysub);
        bt_interp = single(interp2(single(b_xgrid),single(b_ygrid),double(bed_sourcefilt),melmask.dated(p).x,melmask.dated(p).y,'cubic'));
        [xt,yt] = polyxpoly(S.X,S.Y,melmask.dated(p).x,melmask.dated(p).y); %intersect the centerline & terminus & grab depth there
        if sum(~isnan(bt_interp)) > 0
            Ht_float = -(rho_w/rho_i)*bt_interp;
            Hc_float = -(rho_w/rho_i)*single(interp2(single(b_xgrid),single(b_ygrid),double(bed_sourcefilt),xt(1),yt(1),'nearest'));
        else
            bt_interp = single(interp2(single(b_xgrid),single(b_ygrid),double(bed_errfilt),melmask.dated(p).x,melmask.dated(p).y,'cubic'));
            Ht_float = -(rho_w/rho_i)*bt_interp;
            Hc_float = -(rho_w/rho_i)*single(interp2(single(b_xgrid),single(b_ygrid),double(bed_sourcefilt),xt(1),yt(1),'nearest'));
        end

        %calculate the ice thickness for icebergs under the assumption of flotation
        floating = rho_w/(rho_w-rho_i); %multiplier to convert elevation to thickness
        Hi_float = floating.*melange;

        % %plot the DEM if you want to check data coverage
        % figure; set(gcf,'position',[50 50 800 800]);
        % imagesc(M.DEM.x,M.DEM.y,melange); axis xy equal; hold on;
        % melange_cmap = cmocean('thermal',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap);
        % set(gca,'clim',[0 16],'fontsize',16); cbar = colorbar('fontsize',16); cbar.Label.String = 'elevation (m a.s.l.)';
        % set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]);
        % xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
        % set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',16);
        % xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);

        %filter out negative thicknesses
        Hi_float(Hi_float<0) = NaN; Ht_float(Ht_float<0) = NaN; Hc_float(Hc_float<0) = NaN; 

        %only save data from the "terminus" if the inland boundary is
        %actually the terminus and not just the DEM edge
        if MP(i).T.qualflag(p) == 1
            Hmax = [Hmax; max(Hi_float,[],"all"), prctile(Ht_float,[95],"all"), Hc_float];
            Hiqr = [Hiqr; prctile(Ht_float,[25 75],"all")'];
        else
            Hmax = [Hmax; max(Hi_float,[],"all"), NaN, NaN];
            Hiqr = [Hiqr; NaN, NaN];
        end
        
        clear M melange H*_float BB *_xgrid *_ygrid *_interp xt yt;
    end
    disp(['Maximum iceberg thickness range: [',num2str(round(min(Hmax(:,2:3),[],"all"))),...
            '-',num2str(round(max(Hmax(:,2:3),[],"all"))),' m]']);
    disp(['Iceberg common thickness: [',num2str(round(nanmedian(Hiqr(:,1)))),...
            '-',num2str(round(nanmedian(Hiqr(:,2)))),' m]']);
    H_maxes(i,:) = [max(Hmax(:,1)),max(Hmax(:,2)),max(Hmax(:,3))];
    H_range(i,:) = [nanmedian(Hiqr(:,1)),nanmedian(Hiqr(:,2)),max(Hmax(:,3))];

    clear site_abbrev mask_dates S BB bed_* DEM_*
    cd ../..
    %save data to a temporary file in case it crashes
    save('Greenland_iceberg_buttressing.mat',"H_range",'H_maxes','-mat');
end

%% Calculate buttressing 

disp('Buttressing force (x10^7 N/m) for...')
for j = 1:length(H_range)
    site_abbrev = char(folderNames(i));
    for k = length(ang_range)
        for l = 1:3
            [taus, gammas] = calculate_taus(H_range(j,l),ang_range(k),E,rho_i,rho_w);
            Fbutt_bottom(j,l) = -sum(taus)/(gammas(3)*cosd(ang_range(k)));
            Fbutt_top(j,l) = -sum(taus)/((H_range(j,l)*(cosd(ang_range(k)) - E*sind(ang_range(k)))) - gammas(3)*cosd(ang_range(k)));
        end
    end
    disp([site_abbrev,': H=',num2str(H_range(j,l)),' & ang=',num2str(ang_range(k)),': ',num2str(round(Fbutt_bottom(j,l)/10^7,2))]);
end



function [taus,gammas] = calculate_taus(H,ang,E,rho_i,rho_w)

gammas(2) = NaN;
gammas(3)= H*(1 - rho_i/rho_w - (E/2)*tand(ang));
gammas(1) = gammas(3)+E*H*tand(ang);

taus(1) = (1/12)*rho_w*9.81*cosd(ang)*(H-gammas(1))^2*(H+2*gammas(1));
taus(2) = (1/12)*rho_w*9.81*sind(ang)*E^3*H^3;
taus(3) = -(1/12)*rho_w*9.81*cosd(ang)*(H-gammas(3))^2*(H+2*gammas(3));
end
