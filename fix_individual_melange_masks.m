function [removed_flag] = fix_individual_melange_masks(root_dir,site_abbrev,melmask,p)
%%% Use this code to fix one-off strange melange masks that you may suspect
%%% when analyzing melange distributions or melange velocity coherence
% clearvars; close all;
%
% %make sure the cmocean package and the wgs2ps code are in the path
% addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
%     '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/');
elev_cmap = cmocean('thermal',1001); elev_cmap(1,:) = [1 1 1];
%
% %specify site
% site_abbrev = 'ASG'; %3-letter abbrevation
% root_dir = '/Volumes/Jokulhaup_5T/Greenland-melange/';
% cd([root_dir,site_abbrev,'/']);

% %specify the index for the date of interest
% p=29; %example for ASG 20150903
maskref = p;

%load the data
load([root_dir,site_abbrev,'/',site_abbrev,'-melange-masks.mat']);
load([root_dir,site_abbrev,'/DEMs/',site_abbrev,'-',melmask.dated(p).datestring,'_melange-DEMfilled.mat']);

%plot the filled DEM created with the old mask
figure1 = figure; set(gcf,'position',[50 50 1600 600]);
imagesc(double(M.DEM.x),double(M.DEM.y),double(M.DEM.z_filled)); axis xy equal; hold on;
colormap(gca,elev_cmap); set(gca,'clim',[0 80]); cbar = colorbar;
plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',3); hold on;
title(melmask.dated(p).datestring);

%make sure the DEM covers a good area or if the area is too small anyway
%check that the starting vertex for the melange mask is in the ocean
size_answer = questdlg('Does the DEM cover a few square kilometers?',...
    'DEM size check','Yes','No','Yes');
switch size_answer
    case 'Yes'
        
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
        termzoom_question = questdlg('Zoom in to better trace just seaward of the terminus?',...
            'Terminus Zoom Option','1) Yes!','2) No!','1) Yes!');
        figure(figure1);
        switch termzoom_question
            case '1) Yes!'
                disp('Click on UL & LR corners of a box bounding the terminus to zoom in');
                [a] = ginput(2);
                set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                drawnow; clear a;
            case '2) No!'
                %do nothing!
        end
        clear termzoom_question;
        disp('Trace just seaward along the terminus, intersecting each side of the fjord');
        % [term_x,term_y,~] = improfile; %old approach to trace the terminus
        pline = drawpolyline('Color','c','Linewidth',1);
        term_x = pline.Position(:,1); term_y = pline.Position(:,2); clear pline;
        [Z.term.x,Z.term.y] = poly2cw(term_x,term_y);
        
        %find the intercepts with the fjord mask
        out_intercept = []; out_interceptx = []; out_intercepty = [];
        for i = 1:length(outline_x)-1
            [xi,yi] = polyxpoly(outline_x(i:i+1),outline_y(i:i+1),Z.term.x,Z.term.y); %find the intersections of the terminus trace with the melange outline
            if ~isempty(xi)
                out_intercept = [out_intercept i]; out_interceptx = [out_interceptx xi]; out_intercepty = [out_intercepty yi];
            end
            clear xi yi;
        end
        
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
        disp('Recreating the melange mask... this will take ~1 hr!');
        in = inpolygon(ZXgrid,ZYgrid,melpoly_x,melpoly_y);
        M.mask.fjord = zeros(size(M.DEM.z_filled));
        M.mask.fjord(in) = 1;
        M.mask.fjord = round(M.mask.fjord); M.mask.fjord = logical(M.mask.fjord);
        M.DEM.z_filled(M.mask.fjord~=1) = NaN;
        save([root_dir,site_abbrev,'/DEMs/',site_abbrev,'-',melmask.dated(p).datestring,'_melange-DEMfilled.mat'],'M','m','-v7.3');
        
        %set output
        removed_flag = 'fixed';
        
    case 'No'
        %remove from the melange mask matfile
        mel.uncropped = melmask.uncropped;
        mel.dated = melmask.dated(1:p-1);
        mel.dated(p:length(melmask.dated)-1) = melmask.dated(p+1:length(melmask.dated));
        clear melmask; melmask = mel; clear mel;
        save([root_dir,site_abbrev,'/',site_abbrev,'-melange-masks.mat'],'melmask','-v7.3');
        
        %delete the filled DEM & iceberg distribution files (if they exist)
        recycle('on'); 
        disp('Deleting the DEM and associated size distribution files...');
        delete([root_dir,site_abbrev,'/DEMs/',site_abbrev,'-',melmask.dated(p).datestring,'_melange-DEMfilled.mat']);
        %delete size distribution files if they exist
        if exist([root_dir,site_abbrev,'/',site_abbrev,'-',bad_date,'-iceberg-distribution.csv']) == 2
            delete([root_dir,site_abbrev,'/',site_abbrev,'-',bad_date,'-iceberg-distribution.csv']);
            delete([root_dir,site_abbrev,'/',site_abbrev,'-',bad_date,'-iceberg-distribution-subsets.csv']);
            delete([root_dir,site_abbrev,'/',site_abbrev,'-',bad_date,'_melange-distribution_subplots.png']);
        end
        %delete automatically-fit and manually-adjusted modeled size
        %distribution files if they exist
        if exist([root_dir,site_abbrev,'/models/',site_abbrev,'-',bad_date,'_model.png']) == 2
            delete([root_dir,site_abbrev,'/models/',site_abbrev,'-',bad_date,'_model.png']);
            disp('Delete the corresponding row in the time-stamped size distribution model parameters CSV and resave');
        end
        if exist([root_dir,site_abbrev,'/manually_adjusted_models/',site_abbrev,'-',bad_date,'-parameters-adjusted.csv']) == 2
            delete([root_dir,site_abbrev,'/manually_adjusted_models/',site_abbrev,'-',bad_date,'-parameters-adjusted.csv']);
            delete([root_dir,site_abbrev,'/manually_adjusted_models/',site_abbrev,'-',bad_date,'-model-adjusted.png']);
        end
        
        %prompt the user to remove the dated elevations & derived data
        %products if they exist
        disp('If you have already extracted size distributions & elevation profiles, delete the appropriate column from:');
        disp([site_abbrev,'_centerline_elevations.csv']);
        disp([site_abbrev,'_transect_elevations.csv']);
        
        %set output
        removed_flag = 'removed';
end

end