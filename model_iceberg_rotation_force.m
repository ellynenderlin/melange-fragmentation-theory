% calculate theoretical buttressing needed for icebergs trying to calve
% from termini & with varying thickness and rotation angles (keep aspect
% ratio fixed for simplicity)
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/');
addpath('/Users/ellynenderlin/Research/miscellaneous/melange-fragmentation-code/');

%specify variable ranges, paths, and constants
output_dir ='/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';
% H_range = [100:50:1000]; %range of iceberg thicknesses
ang_range = [0.1,1]; %range of angles from vertical
Ef = 0.31; %width-to-thickness ratio for floating terminus
Eg = 0.31; %width-to-thickness ratio for grounded terminus
rho_i = 900; rho_w = 1026; %densities


%% ONLY RUN ONCE: Estimate thicknesses using bed and surface elevations
H_maxes = []; %dummy variables to hold maximum iceberg thickness estimates
H_range = []; %dummy variables to hold quartiles & max terminus thickness estimates

% read in the melange data
load([output_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']);
for i = 1:length(MP)
    mel(i).termflag = MP(i).Z.termflag;
end
clear MP;

% loop through each site and grab floatation-based thickness estimates
cd(output_dir);
contents = dir(output_dir);
folders = contents([contents.isdir]);
validFolders = folders(~ismember({folders.name}, {'.', '..'}));
folderNames = {validFolders.name};
% load('GrIS-iceberg_buttressing.mat'); %reload data

for i = 1:length(folderNames)
    site_abbrev = char(folderNames(i));
    cd([output_dir,site_abbrev]);
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
        if mel(i).termflag(p) == 1
            mask_xext(p,:) = [min(melmask.dated(p).x), max(melmask.dated(p).x)];
            mask_yext(p,:) = [min(melmask.dated(p).y), max(melmask.dated(p).y)];
        else
            mask_xext(p,:) = [NaN, NaN];
            mask_yext(p,:) = [NaN, NaN];
        end
    end
    cd('shapefiles/')
    S = shaperead([site_abbrev,'-centerline.shp']);

    %read in BedMachine Greenland v5 (elevations referenced to the geoid)
    cd('/Users/ellynenderlin/Research/miscellaneous/elevations/BedMachine/');
    bx = ncread('BedMachineGreenland-v5.nc','x');
    by = ncread('BedMachineGreenland-v5.nc','y');
    bz = ncread('BedMachineGreenland-v5.nc','bed'); %bed elevations
    bze = ncread('BedMachineGreenland-v5.nc','errbed'); %bed elevations
    bzs = ncread('BedMachineGreenland-v5.nc','source'); %data interp (0-2 = direct & mass conservation)
    % bzi = ncread('BedMachineGreenland-v5.nc','dataid'); %data source (>2 = seismic or multibeam bathy)
    %bg = ncread('BedMachineGreenland-v5.nc','geoid'); %geoid elevations relative to WGS84 ellipsoid
    cd([output_dir,site_abbrev]);

    %subset the bedmap to the melange area
    % BB = [min(melmask.uncropped.x),min(melmask.uncropped.y); min(melmask.uncropped.x),max(melmask.uncropped.y);...
    %     max(melmask.uncropped.x),max(melmask.uncropped.y); max(melmask.uncropped.x),min(melmask.uncropped.y)];
    BB = [min(mask_xext(:,1))-200,min(mask_yext(:,1)-200); min(mask_xext(:,1))-200,max(mask_yext(:,2)+200);...
        max(mask_xext(:,2)+200),max(mask_yext(:,2)+200); max(mask_xext(:,2)+200),min(mask_yext(:,1)-200)];
    bed_xind = [find(bx<=min(BB(:,1)),1,'last'),find(bx>=max(BB(:,1)),1,'first')];
    bed_yind = [find(by>=max(BB(:,2)),1,'last'),find(by<=min(BB(:,2)),1,'first')];
    bed_xsub = bx(min(bed_xind):max(bed_xind)); bed_ysub = by(min(bed_yind):max(bed_yind));
    bed_zsub = bz(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))'; %rotate so rows are y, columns are x
    bed_zesub = bze(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))';
    bed_sourcesub = bzs(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))';
    % bed_idsub = bzi(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))';
    %bed_gsub = bg(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))'; %rotate so rows are y, columns are x
    
    %filter wonky data from bad interpolation & plot
    bed_sourcefilt = bed_zsub; bed_sourcefilt(bed_sourcesub>2) = NaN;
    bed_errfilt = bed_zsub; bed_errfilt(bed_zesub>=200) = NaN;
    [b_xgrid,b_ygrid] = meshgrid(bed_xsub,bed_ysub);

    %plot to check it looks reasonable
    % figure; imagesc(bed_xsub,bed_ysub,bed_zsub); axis xy equal; colormap jet; colorbar; hold on; %plot to check subsetting
    % for j = 1:length(melmask.dated)
    %     plot(melmask.dated(j).x,melmask.dated(j).y,'k'); hold on;
    % end
    % drawnow;
    clear bg bx by bz*;

    %loop through the dates to grab terminus thickness estimates
    cd([output_dir,site_abbrev,'/DEMs/']);
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
        clear M;

        %get rid of the melange polygon points that are not over the terminus
        [~,on] = inpolygon(melmask.dated(p).x,melmask.dated(p).y,melmask.uncropped.x,melmask.uncropped.y);
        melmask.dated(p).x(on) = []; melmask.dated(p).y(on) = [];
        clear on;

        %interpolate to the same grid
        % [z_xgrid,z_ygrid] = meshgrid(M.DEM.x,M.DEM.y);
        bt_interp = single(interp2(single(b_xgrid),single(b_ygrid),double(bed_sourcefilt),melmask.dated(p).x,melmask.dated(p).y,'linear'));
        [xt,yt] = polyxpoly(S.X,S.Y,melmask.dated(p).x,melmask.dated(p).y); %intersect the centerline & terminus & grab depth there
        if sum(~isnan(bt_interp)) > 0
            Ht_float = -(rho_w/rho_i)*bt_interp;
            if ~isempty(xt)
                Hc_float = -(rho_w/rho_i)*single(interp2(single(b_xgrid),single(b_ygrid),double(bed_sourcefilt),xt(1),yt(1),'nearest'));
            else
                Hc_float = NaN;
            end
        else
            bt_interp = single(interp2(single(b_xgrid),single(b_ygrid),double(bed_errfilt),melmask.dated(p).x,melmask.dated(p).y,'linear'));
            Ht_float = -(rho_w/rho_i)*bt_interp;
            if ~isempty(xt)
                Hc_float = -(rho_w/rho_i)*single(interp2(single(b_xgrid),single(b_ygrid),double(bed_errfilt),xt(1),yt(1),'nearest'));
            else
                Hc_float = NaN;
            end
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
        if mel(i).termflag(p) == 1
            Hmax = [Hmax; max(Hi_float,[],"all"), prctile(Ht_float,[95],"all"), Hc_float];
            Hiqr = [Hiqr; prctile(Ht_float,[25 75],"all")'];
        else
            Hmax = [Hmax; max(Hi_float,[],"all"), NaN, NaN];
            Hiqr = [Hiqr; NaN, NaN];
        end
        save([output_dir,site_abbrev,'/',site_abbrev,'-iceberg_thicknesses.mat'],"Hmax",'Hiqr','-mat');
        
        % clear z_xgrid z_ygrid;
        clear melange H*_float BB *_interp xt yt;
    end
    disp(['Maximum iceberg thickness range: [',num2str(round(min(Hmax(:,2:3),[],"all"))),...
            '-',num2str(round(max(Hmax(:,2:3),[],"all"))),' m]']);
    disp(['Iceberg common thickness: [',num2str(round(nanmedian(Hiqr(:,1)))),...
            '-',num2str(round(nanmedian(Hiqr(:,2)))),' m]']);
    H_maxes(i,:) = [max(Hmax(:,1)),max(Hmax(:,2)),max(Hmax(:,3))];
    H_range(i,:) = [nanmedian(Hiqr(:,1)),nanmedian(Hiqr(:,2)),max(Hmax(:,3))];

    clear site_abbrev mask_dates S BB bed_* DEM_* *_xgrid *_ygrid mask_*ext;
    cd ../..
    %save data to a file
    save('GrIS-iceberg_buttressing.mat',"H_range",'H_maxes','-mat');
end
clear mel;

%% Calculate buttressing 
warning off; close all; drawnow;

%load the data
cd(output_dir);
contents = dir(output_dir);
folders = contents([contents.isdir]);
validFolders = folders(~ismember({folders.name}, {'.', '..'}));
folderNames = {validFolders.name};
load([output_dir,'GrIS-melange_centerline-elev-speed-terminus.mat']);
load([output_dir,'GrIS-iceberg_buttressing.mat']);

%define seasons
seasons = [12,1,2;3,4,5;6,7,8;9,10,11]; season_names = {'DJF','MAM','JJA','SON'};

%define the fraction of the ice thickness relative to just barely grounded
float_frac = [0.5,0.9,1];

% %use 95th percentile thickness along terminus since some
% %glaciers don't have a centerline max thickness value
% H_range(:,3) = H_maxes(:,2);

%initialize a figure
figure; 
mo_cmap = cmocean('phase',12); close all;

disp('Buttressing force (x10^6 N/m) for...')
for i = 1:length(MP)
    site_abbrev = char(folderNames(i));
    load([output_dir,site_abbrev,'/',site_abbrev,'-iceberg_thicknesses.mat'])

    %isolate the dates for data sorting
    for p = 1:length(MP(i).Z.date)
        term_trace(p) = MP(i).Z.termflag(p);
        % zdate(p) = convert_to_decimaldate(char(MP(i).Z.date(p)));
        datest = datetime(MP(i).Z.date{p},'InputFormat','yyyyMMdd');
        zyrs(p) = year(datest); zmos(p) = month(datest);
        clear datest;
    end

    %create dummy vectors to hold variables for the CSV output
    min_ang = NaN(4,9); max_ang = NaN(4,9); 

    %group by season
    for p = 1:4
        seas_inds = find(ismember(zmos,seasons(p,:))==1);
        if ~isempty(seas_inds)
            H_range(p,:) = [nanmean(Hiqr(seas_inds,1)), nanmean(Hiqr(seas_inds,2)), nanmean(Hmax(seas_inds,2))]; %25th, 75th, 95th percentiles
        else
            H_range(p,:) = [NaN, NaN, NaN]; %no data for that season
        end
        clear seas_inds;

        %loop through the angles of consideration
        for j = 1:2
            disp(['ANGLE = ',num2str(ang_range(j))])

            %loop through the flotation fractions
            for k = 1:length(float_frac)
                if float_frac(k) < 1
                    calve_mode = 'floating';
                    E = Ef;
                else
                    calve_mode = 'grounded';
                    E = Eg;
                end

                %calculate buttressing needed to resist rotation of icebergs of
                %different sizes
                for l = 1:3
                    Fbutt(k,l) = calculate_taus(float_frac(k)*H_range(p,l),ang_range(j),E,rho_i,rho_w,calve_mode);
                end
                %display the maximum buttressing needed based on floatation &
                %bottom-out rotation
                disp([site_abbrev,': floatation fraction = ',num2str(float_frac(k))]);
                disp(['maximum bottom-out F_{buttressing}: ',num2str(round(Fbutt(k,l)/10^6,2))])
                %display typical buttressing IQR for the terminus based on floatation
                disp([' typical bottom-out F_{buttressing} IQR: ',num2str(round(Fbutt(k,1)/10^6,2)),'-',num2str(round(Fbutt(k,2)/10^6,2))]);
                % disp([' typical top-out F_{buttressing} IQR: ',num2str(round(Fbutt_top(k,1)/10^6,2)),'-',num2str(round(Fbutt_top(k,2)/10^6,2))]);

            end

            %add data to the matrices used to make CSVs
            if  j == 1 %smaller angle
                for k = 1:length(float_frac)
                    min_ang(p,3*(k-1)+1:3*(k-1)+3) = round(Fbutt(k,:)/10^6,2);
                    if k == 1
                        Fbutt_thinminf(i,p,1:3) = min_ang(p,3*(k-1)+1:3*(k-1)+3);
                    elseif k == 2
                        Fbutt_minf(i,p,1:3) = min_ang(p,3*(k-1)+1:3*(k-1)+3);
                    elseif k == 3
                        Fbutt_ming(i,p,1:3) = min_ang(p,3*(k-1)+1:3*(k-1)+3);
                    end
                    % %plot the data
                    % if k == 1
                    %     % errorbar(i,nanmean(min_ang(p,3*(k-1)+1:3*(k-1)+2)),min_ang(p,3*(k-1)+1),min_ang(p,3*(k-1)+2),...
                    %     %     'o','color',mo_cmap(p*3-2,:),'linewidth',1.5); hold on;
                    % elseif k == 2
                    %     errorbar(i,nanmean(min_ang(p,3*(k-1)+1:3*(k-1)+2)),min_ang(p,3*(k-1)+1),min_ang(p,3*(k-1)+2),...
                    %         's','color',mo_cmap(p*3-2,:),'linewidth',1.5); hold on;
                    % else
                    %     errorbar(i,nanmean(min_ang(p,3*(k-1)+1:3*(k-1)+2)),min_ang(p,3*(k-1)+1),min_ang(p,3*(k-1)+2),...
                    %         'd','color',mo_cmap(p*3-2,:),'linewidth',1.5); hold on;
                    % end
                end
            else %larger angle
                for k = 1:length(float_frac)
                    max_ang(p,3*(k-1)+1:3*(k-1)+3) = round(Fbutt(k,:)/10^6,2);
                    if k == 1
                        Fbutt_thinmaxf(i,p,1:3) = max_ang(p,3*(k-1)+1:3*(k-1)+3);
                    elseif k == 2
                        Fbutt_maxf(i,p,1:3) = max_ang(p,3*(k-1)+1:3*(k-1)+3);
                    elseif k == 3
                        Fbutt_maxg(i,p,1:3) = max_ang(p,3*(k-1)+1:3*(k-1)+3);
                    end
                end
            end

            clear Fbutt;
        end

        %average the seasonal melange buttressing estimates across all years
        Fb_Meng(p,1) = round(nanmean(MP(i).B.butt_Meng(4,p,:))/10^6,2); %4th row is for icebergs with freeboard > 3m
        Fb_Amundson(p,1) = round(nanmean(MP(i).B.butt_Amundson(4,p,:))/10^6,2); %4th row is for icebergs with freeboard > 3m
        Fbutt_obs(i,p,:) = [Fb_Meng(p,1), Fb_Amundson(p,1)];

        % %add the estimated buttressing to the plot
        % plot(i,Fb_Meng(p,1),'+','color',mo_cmap(p*3-2,:),'linewidth',2,'markersize',16); hold on;
        % plot(i,Fb_Amundson(p,1),'x','color',mo_cmap(p*3-2,:),'linewidth',2,'markersize',16); hold on;
        % drawnow;
    end
    clear H_range Hmax Hiqr;


    %combine the results for the different angles & the estimates of actual
    %melange buttressing in a single CSV
    Tnames = table(char(season_names')); Tnames.Properties.VariableNames = "Season"; 
    T = array2table([Fb_Meng,Fb_Amundson,min_ang,max_ang]);
    column_names = ["Buttressing F/W - Meng (MN/m)","Buttressing F/W - Amundson (MN/m)",...
        "Torque (MN/m; 0.1 deg rot, 50% of 25th H)", "Torque (MN/m; 0.1 deg rot, 50% of 75th H)", "Torque (MN/m; 0.1 deg rot, 50% of 95th H)",...
        "Torque (MN/m; 0.1 deg rot, 90% of 25th H)", "Torque (MN/m; 0.1 deg rot, 90% of 75th H)", "Torque (MN/m; 0.1 deg rot, 90% of 95th H)",...
        "Torque (MN/m; 0.1 deg rot, 100% of 25th H)", "Torque (MN/m; 0.1 deg rot, 100% of 75th H)", "Torque (MN/m; 0.1 deg rot, 100% of 95th H)",...
        "Torque (MN/m; 1.0 deg rot, 50% of 25th H)", "Torque (MN/m; 1.0 deg rot, 50% of 75th H)", "Torque (MN/m; 1.0 deg rot, 50% of 95th H)",...
        "Torque (MN/m; 1.0 deg rot, 90% of 25th H)", "Torque (MN/m; 1.0 deg rot, 90% of 75th H)", "Torque (MN/m; 1.0 deg rot, 90% of 95th H)",...
        "Torque (MN/m; 1.0 deg rot, 100% of 25th H)", "Torque (MN/m; 1.0 deg rot, 100% of 75th H)", "Torque (MN/m; 1.0 deg rot, 100% of 95th H)"];
    T.Properties.VariableNames = column_names; %T.Properties.RowNames = season_names';
    Tfinal = [Tnames,T];
    writetable(Tfinal,[output_dir,site_abbrev,'/',site_abbrev,'-buttressing-forces.csv']);
    disp([site_abbrev, ' buttressing text file written']);
    disp(' ')

    clear Tnames T Tfinal min_ang max_ang Fb_* term_trace zyrs zmos;
end
%format the data plot
set(gcf,'position',[50 50 1200 1200]);
sub1 = subplot(3,1,1); 
fill([0,5,0,0],[0,5,5,0],'k','FaceAlpha',0.2,'EdgeColor','none'); hold on;
plot([0,5],[0,5],'--k','linewidth',1.5); hold on;
text(1.05,1.2,'unbuttressed','fontsize',16,'rotation',22.5); 
text(1.15,1.0,'buttressed','fontsize',16,'rotation',22.5); 
sub2 = subplot(3,1,2);
fill([0,5,0,0],[0,5,5,0],'k','FaceAlpha',0.2,'EdgeColor','none'); hold on;
text(1.05,1.2,'unbuttressed','fontsize',16,'rotation',22.5); 
text(1.15,1.0,'buttressed','fontsize',16,'rotation',22.5);  
plot([0,5],[0,5],'--k','linewidth',1.5); hold on;
sub3 = subplot(3,1,3);
fill([0,5,0,0],[0,5,5,0],'k','FaceAlpha',0.2,'EdgeColor','none'); hold on;
text(1.05,1.2,'unbuttressed','fontsize',16,'rotation',22.5); 
text(1.15,1.0,'buttressed','fontsize',16,'rotation',22.5);    
plot([0,5],[0,5],'--k','linewidth',1.5); hold on;
for p = 1:4
    ymid = nanmean(Fbutt_obs(:,p,:),3); yneg_err = ymid - min(Fbutt_obs(:,p,:),[],3); ypos_err = min(Fbutt_obs(:,p,:),[],3) - ymid; 
    
    %needed vs observed: barely grounded (GL)
    subplot(sub1);
    xmid = nanmean(Fbutt_ming(:,p,1:2),3); xneg_err = xmid - min(Fbutt_ming(:,p,1:2),[],3); xpos_err = min(Fbutt_ming(:,p,1:2),[],3) - xmid; 
    pl(p) = errorbar(ymid,xmid,xneg_err,xpos_err,yneg_err,yneg_err,'s','color',mo_cmap(p*3-2,:),'linewidth',2); hold on;
    plot(ymid,Fbutt_ming(:,p,3),'*','color',mo_cmap(p*3-2,:),'linewidth',2); hold on;
    clear xmid;
    % xmid = nanmean(Fbutt_maxg(:,p,1:2),3); xneg_err = xmid - min(Fbutt_maxg(:,p,1:2),[],3); xpos_err = min(Fbutt_maxg(:,p,1:2),[],3) - xmid; 
    % errorbar(xmid,ymid,yneg_err,yneg_err,xneg_err,xpos_err,'d','color',mo_cmap(p*3-2,:),'linewidth',1); hold on;
    % plot(Fbutt_maxg(:,p,3),ymid,'+','color',mo_cmap(p*3-2,:),'linewidth',1); hold on;
    % clear xmid;

    %needed vs observed: 90% GL thickness
    subplot(sub2);
    xmid = nanmean(Fbutt_minf(:,p,1:2),3); xneg_err = xmid - min(Fbutt_minf(:,p,1:2),[],3); xpos_err = min(Fbutt_minf(:,p,1:2),[],3) - xmid; 
    errorbar(ymid,xmid,xneg_err,xpos_err,yneg_err,yneg_err,'s','color',mo_cmap(p*3-2,:),'linewidth',2); hold on;
    plot(ymid,Fbutt_minf(:,p,3),'*','color',mo_cmap(p*3-2,:),'linewidth',2); hold on;
    clear xmid;
    % xmid = nanmean(Fbutt_maxf(:,p,1:2),3); xneg_err = xmid - max(Fbutt_minf(:,p,1:2),[],3); xpos_err = min(Fbutt_maxf(:,p,1:2),[],3) - xmid; 
    % errorbar(xmid,ymid,yneg_err,yneg_err,xneg_err,xpos_err,'d','color',mo_cmap(p*3-2,:),'linewidth',1); hold on;
    % plot(Fbutt_maxf(:,p,3),ymid,'+','color',mo_cmap(p*3-2,:),'linewidth',1); hold on;
    % clear xmid;
    
    %needed vs observed: 50% GL thickness
    subplot(sub3);
    xmid = nanmean(Fbutt_thinminf(:,p,1:2),3); xneg_err = xmid - min(Fbutt_thinminf(:,p,1:2),[],3); xpos_err = min(Fbutt_thinminf(:,p,1:2),[],3) - xmid; 
    errorbar(ymid,xmid,xneg_err,xpos_err,yneg_err,yneg_err,'s','color',mo_cmap(p*3-2,:),'linewidth',2); hold on;
    plot(ymid,Fbutt_thinminf(:,p,3),'*','color',mo_cmap(p*3-2,:),'linewidth',2); hold on;
    clear xmid;
    
    % %observed for each site
    % subplot(sub1);
    % plot([1:1:length(MP)],Fbutt_obs(:,p,1),'+','color',mo_cmap(p*3-2,:),'linewidth',2,'markersize',16); hold on;
    % plot([1:1:length(MP)],Fbutt_obs(:,p,2),'x','color',mo_cmap(p*3-2,:),'linewidth',2,'markersize',16); hold on;
    % subplot(sub2);
    % plot([1:1:length(MP)],Fbutt_obs(:,p,1),'+','color',mo_cmap(p*3-2,:),'linewidth',2,'markersize',16); hold on;
    % plot([1:1:length(MP)],Fbutt_obs(:,p,2),'x','color',mo_cmap(p*3-2,:),'linewidth',2,'markersize',16); hold on;
end
subplot(sub1); set(gca,'ylim',[0 2.5],'xlim',[0 5],'fontsize',20); 
pl_leg = legend(pl,season_names);
text(0.05,2.2,'a) Grounding Line Thickness (GLH)','fontsize',20); 
subplot(sub2); set(gca,'ylim',[0 2.5],'xlim',[0 5],'fontsize',20); 
rectangle('Position',[4.0 1.75 0.95 0.65],'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5);
errorbar(4.15,2.3,0,0,0.1,0.1,'s','color','k','linewidth',2); hold on;
text(4.3,2.3,'M24-A25 span','fontsize',16);
errorbar(4.15,2.1,0.1,0.1,'s','color','k','linewidth',2); hold on;
text(4.3,2.1,'IQR','fontsize',16);
plot(4.15,1.9,'*','color','k','linewidth',2); hold on;
text(4.3,1.9,'95th-percentile','fontsize',16);
text(0.05,2.2,'b) 90% GLH','fontsize',20); 
ylabel('Bottom-out calving torque (x10^6 N/m)','fontsize',20);
subplot(sub3); set(gca,'ylim',[0 2.5],'xlim',[0 5],'fontsize',20); 
text(0.05,2.2,'c) 50% GLH','fontsize',20); 
xlabel('Observed buttressing (x10^6 N/m)','fontsize',20);
% xticklabels = folderNames;
% set(gca,'xtick',[1:1:length(MP)],'xticklabel',xticklabels,'fontsize',20); 
% ylabel('Buttressing F/W (x10^6 N/m)','fontsize',20);
saveas(gcf,[output_dir,'GrIS-buttressing-intercomparison-scatterplots.png'],'png'); %save the plots


function [Ft] = calculate_taus(H,ang,E,rho_i,rho_w,calve_mode)

% friction coefficients for terminus and fjord bottom
muT=0; muN=0;

% calculate gammas (iceberg length above water)
if contains(calve_mode,'float')
    gammas(3)= H*(1 - rho_i/rho_w - (E/2)*tand(ang));
else %grounded
    gammas(3)=H*(1-(rho_i/rho_w)*secd(ang));
end
gammas(2) = NaN;
gammas(1) = gammas(3)+E*H*tand(ang);

%calculate torques
taus(1) = (1/12)*rho_w*9.81*cosd(ang)*(H-gammas(1))^2*(H+2*gammas(1));
taus(2) = (1/12)*rho_w*9.81*sind(ang)*E^3*H^3;
taus(3) = -(1/12)*rho_w*9.81*cosd(ang)*(H-gammas(3))^2*(H+2*gammas(3));
tauB=sum(taus);

% force due to gravity
Fg=rho_i*9.81*H^2*E;

% force due to buoyancy
Fb=rho_w*9.81*E*H/2*(2*H-gammas(3)-gammas(1));

% normal force on lowest corner of the iceberg
Fn=Fg-Fb;

% torque from normal force on bottom of iceberg
tauN=Fn*H/2*(E*cosd(ang)-sind(ang));

%solve separate conditions for floating vs lightly grounded termini
if contains(calve_mode,'float')
    Ft = -tauB/(gammas(3)*cosd(ang)); %bottom-out rotation
    % Ft = -tauB/((H*(cosd(ang) - E*sind(ang))) - gammas(3)*cosd(ang)); %top-out rotation
else
    Ft=(-tauB-tauN-muN*Fn*H/2*(cosd(ang)+E*sind(ang)))/ ... 
        (H*(cosd(ang)-(rho_i/rho_w)-muT/2*(sind(ang)+E*cosd(ang))));
end


end
