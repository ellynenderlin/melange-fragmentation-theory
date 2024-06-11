%%% Use this code to fix one-off strange melange masks that you may suspect
%%% when analyzing melange distributions or melange velocity coherence
clearvars; close all;

%make sure the cmocean package and the wgs2ps code are in the path
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/');
elev_cmap = cmocean('thermal',1001); elev_cmap(1,:) = [1 1 1];

%specify site
site_abbrev = 'ASG'; %3-letter abbrevation
root_dir = '/Volumes/Jokulhaup_5T/Greenland-melange/';
cd([root_dir,site_abbrev,'/']);

%specify the index for the date of interest
p=29; %example for ASG 20150903
maskref = p;

%load the data
load([root_dir,site_abbrev,'/',site_abbrev,'-melange-masks.mat']);
load([root_dir,site_abbrev,'/DEMs/',site_abbrev,'-',melmask.dated(p).datstring,'_melange-DEMfilled.mat']);

%plot the filled DEM created with the old mask
figure1 = figure; set(gcf,'position',[50 50 1600 600]);
imagesc(double(M.DEM.x),double(M.DEM.y),double(M.DEM.z_filled)); axis xy equal; hold on;
colormap(gca,elev_cmap); set(gca,'clim',[0 80]); cbar = colorbar;
plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
title(DEM_name(1:11));

%make sure the starting vertex for the melange mask is placed out in the ocean
[meledge_x, meledge_y] = poly2cw(melmask.uncropped.x,melmask.uncropped.y); %make sure melange outline is a clockwise polygon
[gris_center_x, gris_center_y] = wgs2ps(-41.2, 76.7); % grab the center of the GrIS in PS coordinates
meledge_dist = sqrt((meledge_x-gris_center_x).^2 + (meledge_y-gris_center_y).^2); %find the distance from the GrIS center to each melange outline vertex
start_vert = find(meledge_dist==max(meledge_dist)); %find the farthest vertex & use that as the start of the outline
plot(meledge_x(start_vert(1)),meledge_y(start_vert(1)),'mx','linewidth',3); drawnow;

%check that the starting vertex for the melange mask is in the ocean
answer = questdlg('Where is the starting vertex for the melange mask (pink X)?',...
    'Mask Start Vertex','Ocean','Glacier','Ocean');
switch answer
    case 'Ocean'
        disp('Melange mask start vertex is correct, carry on!');
    case 'Glacier'
        start_dist = sqrt((meledge_x-meledge_x(start_vert(1))).^2 + (meledge_y-meledge_y(start_vert(1))).^2);
        clear start_vert;
        start_vert = find(start_dist==max(start_dist)); %find the farthest vertex from the incorrect default starting vertex
        plot(meledge_x(start_vert(1)),meledge_y(start_vert(1)),'c+','linewidth',3); drawnow;
end

%sort melange mask vertices so that they start at the defined starting vertex in the ocean
if start_vert(1) > 1
    outline_x = [meledge_x(start_vert(1):end); meledge_x(1:start_vert(1))];
    outline_y = [meledge_y(start_vert(1):end); meledge_y(1:start_vert(1))];
else
    outline_x = meledge_x; outline_y = meledge_y;
end

%clear the existing dated melange mask & retrace the terminus
melmask.dated(maskref).x = []; melmask.dated(maskref).y = [];
[term_x,term_y,~] = improfile;
[Z.term.x,Z.term.y] = poly2cw(term_x,term_y);

%recrop the melange mask
term_in = inpolygon(Z.term.x,Z.term.y,outline_x,outline_y); termx = Z.term.x(term_in); termy = Z.term.y(term_in);
termdist = sqrt((outline_x(out_intercept(1))-termx).^2 + (outline_y(out_intercept(1))-termy).^2);
term_start = find(termdist == min(termdist)); term_end = find(termdist == max(termdist));
if term_start > term_end %start at the end of the terminus trace and move to the beginning
    melpoly_x = [outline_x(1:out_intercept(1)); out_interceptx(1); termx(end:-1:1); out_interceptx(end); outline_x(out_intercept(end)+1:end)];
    melpoly_y = [outline_y(1:out_intercept(1)); out_intercepty(1); termy(end:-1:1); out_intercepty(end); outline_y(out_intercept(end)+1:end)];
else %start at the beginning of the terminus trace and move to the end
    melpoly_x = [outline_x(1:out_intercept(1)); out_interceptx(1); termx(1:1:end); out_interceptx(end); outline_x(out_intercept(end)+1:end)];
    melpoly_y = [outline_y(1:out_intercept(1)); out_intercepty(1); termy(1:1:end); out_intercepty(end); outline_y(out_intercept(end)+1:end)];
end

%remake the fjord mask & save the dated melange outline and fjord mask
plot(melpoly_x,melpoly_y,'--c','linewidth',2); hold on;
melmask.dated(maskref).x = melpoly_x; melmask.dated(maskref).y = melpoly_y;
save([root_dir,site_abbrev,'/',site_abbrev,'-melange-masks.mat'],'melmask','-v7.3');
[ZXgrid,ZYgrid] = meshgrid(M.DEM.x,M.DEM.y);
in = inpolygon(ZXgrid,ZYgrid,melpoly_x,melpoly_y);
M.mask.fjord = zeros(size(M.DEM.z_filled));
M.mask.fjord(in) = 1;
M.mask.fjord = round(M.mask.fjord); M.mask.fjord = logical(M.mask.fjord);
save([root_dir,site_abbrev,'/DEMs/',site_abbrev,'-',melmask.dated(p).datestring,'_melange-DEMfilled.mat'],'M','m','-v7.3');


