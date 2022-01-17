%plot AGU 2021 fragmentation figures
clear all; close all;
cd /Users/icebergs/iceberg-fragmentation
site_names = ['HM';'KO';'AG';'IG';'UN';'US';'IB';'UM';'RI';'JI';'KB';'HH';'MG';'KL';'MD';'DJ';'ZI'];
reg_flags = [3;3;3;3;3;3;2;2;2;2;1;1;1;1;4;4;5];
reg_colors = [215,25,28; 253,174,97; 255,255,191; 171,217,233; 44,123,182]/255;
% addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean');
addpath('/Users/icebergs/general-code/');
addpath('/Users/icebergs/general-code/cmocean/');
years = [2011:1:2020];
cmap = cmocean('haline',length(years));
rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)
WV_sigma = 2.9; %DEM uncertainty in meters (used when filtering bad data)
Wmax = 1000; %maximum likely iceberg width (based on qualitative inspection of images during terminus mapping)
Hmax = 800; %threshold thickness in meters that you do not expect icebergs to exceed (grounding line thickness is a good proxy)


%% plot each Landsat images with an overlay of the fjord polygon
for i = 1:length(site_names)
    disp(site_names(i,:));
    cd_to_site = ['cd ',site_names(i,:)]; eval(cd_to_site);
    
    %load & plot landsat panchromatic polar stereo image
    l8 = dir('LC08*'); cd_to_l8 = ['cd ',l8(1).name]; eval(cd_to_l8);
    l8 = dir('LC08*B8PS.TIF');
    [I,S] = geotiffread(l8(1).name);
    x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
    y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
    figure; set(gcf,'position',[50 50 600 600]);
    imagesc(x,y,I); colormap gray; axis xy equal; hold on;
    cd ..
    
    %load and plot the fjord polygon
    melfile = dir('*melange-masks.mat');
    load_mask = ['load ',melfile(1).name]; eval(load_mask);
    %plot dummy lines for colorbar to ensure all years are included
    for j = 1:length(cmap)
       pl(j) = plot(melmask.dated(1).x,melmask.dated(1).y,'-','color',cmap(j,:),'linewidth',1.5); hold on; 
    end
    %plot actual data
    for j = 1:length(melmask.dated)
       plot(melmask.dated(j).x,melmask.dated(j).y,'-','color',cmap(str2num(melmask.dated(j).datestring(1:4))-2011+1,:),'linewidth',1.5); hold on; 
    end
    plot(melmask.uncropped.x,melmask.uncropped.y,'-','color','w','linewidth',2.5); hold on;
    plot(melmask.uncropped.x,melmask.uncropped.y,'-','color',reg_colors(reg_flags(i),:),'linewidth',2); hold on;
    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],...
        'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)],...
        'fontsize',20); grid on;
    xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
    set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000);
    xlabel('Easting (km)','fontsize',20); ylabel('Northing (km)','fontsize',20);
    if reg_flags(i) <=2
        leg = legend(pl,num2str([2011:1:2020]'),'location','west','orientation','vertical');
    else
        leg = legend(pl,num2str([2011:1:2020]'),'location','east','orientation','vertical');
    end
    drawnow;
    saveas(gcf,[site_names(i,:),'_site-map.png'],'png');
    
    cd ..
    disp('moving on...'); close all; drawnow;
end
disp('Finished creating site maps');


%% plot example DEMs
DEMpath = '/Users/icebergs/iceberg-fragmentation/AGU2021/';
maskpath = '/Users/icebergs/iceberg-fragmentation/';

% %uncomment lines below to selectively replot
% clear site_names reg_flags;
% site_names = ['IG';'UN';'IB';'RI';'MD'];
% reg_flags = [3;3;2;2;4];

for i = 1:length(site_names)
    disp(site_names(i,:));
    
    %load the DEM
    cd_to_DEM = ['cd ',DEMpath,site_names(i,:)]; eval(cd_to_DEM);
    DEMs = dir('*dem.tif');
    DEMmat_dates = DEMs(1).name(6:13);
    [Z.z.raw,S] = readgeoraster(DEMs(1).name); %load the DEM
    Z.x = S.XWorldLimits(1)+0.5*S.CellExtentInWorldX:S.CellExtentInWorldX:S.XWorldLimits(2)-0.5*S.CellExtentInWorldX;
    Z.y = S.YWorldLimits(2)-0.5*S.CellExtentInWorldY:-S.CellExtentInWorldY:S.YWorldLimits(1)+0.5*S.CellExtentInWorldY;
    Z.z.raw(Z.z.raw==-9999) = NaN;
    
    %convert to orthometric elevations
    disp('calculating the geoid heights for the DEM...')
    warning off;
    geoid_spacing = 100;
    [ZXgrid,ZYgrid] = meshgrid(Z.x,Z.y);
    for j = 1:geoid_spacing:size(ZXgrid,1)
        for k = 1:geoid_spacing:size(ZXgrid,2)
            [Zlon(ceil(j/geoid_spacing),ceil(k/geoid_spacing)),Zlat(ceil(j/geoid_spacing),ceil(k/geoid_spacing))] = ps2wgs(ZXgrid(j,k),ZYgrid(j,k));
            G(ceil(j/geoid_spacing),ceil(k/geoid_spacing)) = geoidheight(Zlat(ceil(j/geoid_spacing),ceil(k/geoid_spacing)),Zlon(ceil(j/geoid_spacing),ceil(k/geoid_spacing)));
        end
    end
    Z.z.geoid = single(interp2(ZXgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),ZYgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),G,ZXgrid,ZYgrid,'linear'));
    disp('converting from ellipsoidal heights to orthometric elevations...');
    Z.z.ortho = Z.z.raw - Z.z.geoid; %Z.z = rmfield(Z.z,'raw');

    
    %add melange masks
    cd_to_mask = ['cd ',maskpath,site_names(i,:)]; eval(cd_to_mask);
    melfile = dir('*melange-masks.mat');
    load_mask = ['load ',melfile(1).name]; eval(load_mask);
    disp('DEM and melange mask loaded');
    
    %convert melange mask to fjord mask
    in = inpolygon(ZXgrid,ZYgrid,melmask.uncropped.x,melmask.uncropped.y);
    Z.fjord.DEM_maskX = single(melmask.uncropped.x); Z.fjord.DEM_maskY = single(melmask.uncropped.y);
    Z.fjord.DEM_mask = zeros(size(Z.z.ortho));
    Z.fjord.DEM_mask(in) = 1;
    Z.fjord.DEM_mask = round(Z.fjord.DEM_mask); Z.fjord.DEM_mask = logical(Z.fjord.DEM_mask);
    
    %plot DEMs
    [Z,data_mask,gap_mask] = sl_correction(Z,WV_sigma); % perform sea level correction
    M.DEM.x = Z.x; M.DEM.y = Z.y;
    M.DEM.z = Z.z.adjusted; M.DEM.z(Z.z.raw==-9999) = NaN; M.DEM.z(isnan(Z.z.raw)) = NaN; 
    %extract melange elevations
    melange = M.DEM.z;
    melange(isnan(M.DEM.z)) = 0;
    melange(melange<0) = 0; melange(melange>Hmax/(917/(1026-917))) = NaN;
    melange(Z.fjord.DEM_mask==0) = NaN;
    
    %plot the melange DEM
    figDEM = figure; set(figDEM,'position',[550 100 1000 500]);
    imagesc(M.DEM.x,M.DEM.y,melange); axis xy equal;
    melange_cmap = cmocean('thermal',1001); melange_cmap(1,:) = [1 1 1]; colormap(gca,melange_cmap);
    hold on; DEMax = gca;
    set(gca,'clim',[0 16]);  %set(gca,'clim',[0 cmax]);
    cbar = colorbar; cbar.Label.String  = 'elevation (m a.s.l.)';
    if find(nansum(melange)>0,1,'first')-50 > 1; mel_xlims = [find(nansum(melange)>0,1,'first')-50]; else; mel_xlims = 1; end
    if find(nansum(melange)>0,1,'last')+50 < length(M.DEM.x); mel_xlims = [mel_xlims find(nansum(melange)>0,1,'last')+50]; else; mel_xlims = [mel_xlims length(M.DEM.x)]; end
    if find(nansum(melange,2)>0,1,'first')-50 > 1; mel_ylims = [find(nansum(melange,2)>0,1,'first')-50]; else; mel_ylims = 1; end
    if find(nansum(melange,2)>0,1,'last')+50 < length(M.DEM.y); mel_ylims = [mel_ylims find(nansum(melange,2)>0,1,'last')+50]; else; mel_ylims = [mel_ylims length(M.DEM.y)]; end
    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)]); xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
    set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000,'fontsize',28);
    xlabel('Easting (km)','fontsize',28); ylabel('Northing (km)','fontsize',28);
    
    
    %plot dummy lines for colorbar to ensure all years are included
    for j = 1:length(cmap)
        pl(j) = plot(melmask.dated(1).x,melmask.dated(1).y,'-','color',cmap(j,:),'linewidth',1.5); hold on;
    end
    %plot actual data
    for j = 1:length(melmask.dated)
        plot(melmask.dated(j).x,melmask.dated(j).y,'-','color',cmap(str2num(melmask.dated(j).datestring(1:4))-2011+1,:),'linewidth',1.5); hold on;
    end
    plot(melmask.uncropped.x,melmask.uncropped.y,'-','color','w','linewidth',2.5); hold on;
    plot(melmask.uncropped.x,melmask.uncropped.y,'-','color',reg_colors(reg_flags(i),:),'linewidth',2); hold on;
    set(gca,'xlim',[min(melmask.uncropped.x) max(melmask.uncropped.x)],...
        'ylim',[min(melmask.uncropped.y) max(melmask.uncropped.y)],...
        'fontsize',20); grid on;
    xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
    set(gca,'xticklabel',xticks/1000,'yticklabel',yticks/1000);
%     xlabel('Easting (km)','fontsize',28); ylabel('Northing (km)','fontsize',28);
    if reg_flags(i) <=2
        leg = legend(pl,num2str([2011:1:2020]'),'location','northoutside','NumColumns',5);
    else
        leg = legend(pl,num2str([2011:1:2020]'),'location','northoutside','NumColumns',5);
    end
    axpos = get(gca,'position'); set(gca,'position',[axpos(1) 0.11 axpos(3) 0.75]); axpos = get(gca,'position'); 
    legpos = get(leg,'position'); 
    set(leg,'position',[(axpos(1)+axpos(3)/2)-(legpos(3)/2) (axpos(2)+axpos(4))+0.02 legpos(3) legpos(4)]);
    drawnow;
    
    %add title and save
%     title([DEMmat_dates(1:4),'/',DEMmat_dates(5:6),'/',DEMmat_dates(7:8)],'fontsize',20); grid on; drawnow;
    cd_to_AGU = ['cd ',DEMpath]; eval(cd_to_AGU);
    saveas(figDEM,[site_names(i,:),'-',DEMmat_dates(:)','_melange-rawDEM.png'],'png');
    close all; drawnow;
    disp('moving on...');
    clear M Z melange melmask in G Z*grid Zlat Zlon;
end