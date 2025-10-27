function [C,XF,transect_inc] = create_profile_and_transects(varargin)
%SYNTAX (required inputs): [C,XF] = create_profile_and_transects(site_dir,AOI,im_dir,vel_dir)
%SYNTAX (optional inputs must be in this order): [C,XF] = create_profile_and_transects(site_dir,AOI,im_dir,vel_dir,epsg_code,transect_inc,export_type)
%
%DESCRIPTION: Manually delineate a centerline using an area of interest
%(AOI) and an image overlain with a quiver plot of velocities. The AOI,
%image, and velocity data must be in the same projection. Cross-flow
%transects will be produced automatically at a specified uniform increment
%(default = 2000 m), cropped to the extent of the AOI.
%
%INPUT:
%site_dir = root directory for data specific to the study site
%
%AOI = opened shapefile structure, melange mask variable produced by
%create_melange_masks.m (melmask.uncropped), or full path for a shapefile
%that defines the area of interest
%
%im_dir = full path for a Landsat 8 or 9 panchromatic image in the same
%coordinate reference system as the AOI (Landsat default is UTM)
%
%vel_dir = full path for velocity vector geotiffs with 'vx' and 'vy' in
%file names used to identify the two vector components
%
%epsg_code = numeric EPSG code for the coordinate reference system
%(default = 3413 for Greenland Polar Stereographic)
%transect_inc = distance, in meters, between cross-flow transects extending
%from the centerline (default = 2000)
%
%export_type = file type for centerline and transect polylines
%(options: SHP, CSV; default = SHP)
%
%OUTPUT:
%C = X,Y coordinates for the centerline at 2 m increments (adjust 'spacer'
%variable to change the increment) in the same crs as the inputs
%XF = X,Y coordinates for the transects, with a separate structure index
%for each unique transect, in the same crs as the inputs
%SHP or CSV = the centerline and transect polylines will be saved as either
%shapefiles or comma separated value files in a 'shapefiles' directory
%
%DEPENDENCIES:
%interparc.m available through Matlab File Exchange
%(https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc)
%
%EXAMPLE:
%[C,XF] = create_profile_and_transects('/Volumes/Jokulhaup_5T/Greenland-melange/DJG/',melmask,...
%   '/Volumes/Jokulhaup_5T/Greenland-melange/DJG/LC08_L1TP_232009_20130827_20170502_01_T1/',...
%   '/Users/ellynenderlin/Research/miscellaneous/Greenland-VelMosaic_1995-2015/',3413);

%% initialize
close all;
warning('All inputs must be in the same spatial referencing system!');

%required input variables
site_dir = varargin{1}; %output directory
AOI = varargin{2};
im_dir = varargin{3};
vel_dir = varargin{4};

%optional inputs
if nargin >= 5
    epsg_code = varargin{5}; %specified coordinate reference system

    if nargin >= 6
        transect_inc = varargin{6}; %specified transect spacing
        if nargin == 7
            export_type = varargin{7}; %specified file type for polylines
        else
            export_type = 'shp';
        end
    else
        transect_inc = 2000; %default transect spacing (m)
        export_type = 'shp'; %default file type for polyline export
    end
else
    %apply defaults
    epsg_code = 3413; %defaulty crs (Greenland Polar Stereo)
    transect_inc = 2000; %default transect spacing (m)
    export_type = 'shp'; %default file type for polyline export
end
site_crs = projcrs(epsg_code);

spacer = 2; %spacing increment (m): default 2 m is the spatial resolution of WorldView DEMs

%% load data

% read the area of interest shapefile if the AOI isn't provided as a variable
disp('loading the AOI');
if exist('AOI') == 1 %it's a variable in the workspace!
    if isfield(AOI,'X')
        S = AOI;
    elseif isfield(AOI,'uncropped')
        S.X = AOI.uncropped.x; S.Y = AOI.uncropped.y;
    end
elseif exist('AOI') == 7
    if isfile(AOI)
        S = shaperead([site_dir,AOI]);
    elseif ~isfile(AOI) && isfolder([site_dir,'shapefiles'])
        S = shaperead([site_dir,'shapefiles/',AOI]);
    else
        error('cannot locate shapefile in specified directory or shapefiles subdirectory');
    end
else
    error('Unsupported AOI type: Must provide a shapefile path, shapefile variable, or melange mask (melmask.uncropped) variable');
end

%specify Landsat 8/9 image path as an input: use Band 8 for reference image
disp('loading the Landsat panchromatic band');
% im_dir = dir([site_dir,'LC*','/']); %DEFINE AS AN INPUT VARIABLE
ims = dir([im_dir,'L*.TIF']);
for j = 1:length(ims)
    if contains(ims(j).name,'B8')
        ref_image = [ims(j).folder,'/',ims(j).name];
    end
end
clear im_dir ims;
%load the panchromatic Landsat scene & create a glacier/fjord mask
[I,R] = readgeoraster(ref_image);
%             info = geotiffinfo(ref_image);
im.x = R.XWorldLimits(1):R.SampleSpacingInWorldX:R.XWorldLimits(2);
im.y = R.YWorldLimits(2):-R.SampleSpacingInWorldY:R.YWorldLimits(1);
im.z = double(I);
clear I R;

%specify velocity map path as an input
disp('loading the x & y velocity vectors');
% vel_dir = '/Users/ellynenderlin/Research/miscellaneous/Greenland-VelMosaic_1995-2015/'; %DEFINE AS INPUT VARIABLE
vel_tifs = dir([vel_dir,'*.tif']);
%read just the spatial referencing info for the velocity map
for j = 1:length(vel_tifs)
    if contains(vel_tifs(j).name,'vx')
        [~,R] = readgeoraster([vel_tifs(j).folder,'/',vel_tifs(j).name]);
        V.x = single([R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX]);
        V.y = single([R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY]);
        clear R;
    end
end
v_xsub = [find(V.x<=min(S.X),1,'last'),find(V.x>=max(S.X),1,'first')];
v_ysub = [find(V.y>=max(S.Y),1,'last'),find(V.y<=min(S.Y),1,'first')];
%only read a subset of the velocity vector geotiffs
for j = 1:length(vel_tifs)
    if contains(vel_tifs(j).name,'vx')
        I = imread([vel_tifs(j).folder,'/',vel_tifs(j).name],'PixelRegion',{v_ysub,v_xsub});
        %         [I,~] = readgeoraster([vel_tifs(j).folder,'/',vel_tifs(j).name]);
        V.vx = single(I); V.vx(V.vx==-2.0000e+09) = NaN; %velocity in x-direction in m/yr
        clear I R;
    elseif contains(vel_tifs(j).name,'vy')
        I = imread([vel_tifs(j).folder,'/',vel_tifs(j).name],'PixelRegion',{v_ysub,v_xsub});
        %         [I,~] = readgeoraster([vel_tifs(j).folder,'/',vel_tifs(j).name]);
        V.vy = single(I); V.vy(V.vy==-2.0000e+09) = NaN; %velocity in y-direction in m/yr
        clear I;
    end
end
V.x = V.x(min(v_xsub):max(v_xsub)); V.y = V.y(min(v_ysub):max(v_ysub));
[VXgrid,VYgrid] = meshgrid(V.x,V.y);


%% plot the reference image

%crop the image to adjust brightnesses appropriately
xlims = [find(im.x<=min(S.X),1,'last'), find(im.x>=max(S.X),1,'first')];
ylims = [find(im.y>=max(S.Y),1,'last'), find(im.y<=min(S.Y),1,'first')];
im_subset = im.z(min(ylims):max(ylims),min(xlims):max(xlims));
im_subset = im_subset./max(max(im_subset));

%plot the image
fig = figure; set(fig,'position',[50 50 1200 1200]);
imagesc(im.x(min(xlims):max(xlims)),im.y(min(ylims):max(ylims)),imadjust(im_subset)); axis xy equal; colormap gray; drawnow; hold on;

%add the quiver to the map
set(gca,'xlim',[min(S.X) max(S.X)],'ylim',[min(S.Y) max(S.Y)],'fontsize',14);
q = quiver(VXgrid,VYgrid,V.vx,V.vy,'-c'); hold on; %q.AutoScaleFactor = 50;

%overlay the AOI
plot(S.X,S.Y,'-r','linewidth',2); %plot the AOI
drawnow;

%% trace the centerline, fit a spline, and extract points at even increments

%look for an existing centerline file & give the user the option to load it
centerline_file = dir([site_dir,'shapefiles/*-centerline.shp']);
if ~isempty(centerline_file)
    center_reuse = questdlg('Do you want to use the existing centerline (to re-do transect spacing)?',...
        'reuse','1) Yes','2) No','2) No');
    switch center_reuse
        case '1) Yes'
            C = shaperead([centerline_file(1).folder,'/',centerline_file(1).name]);
            figure(fig);
            plot(C.X,C.Y,'-y','linewidth',2); %plot the centerline
            drawnow;
            delineation_flag = 0;
        case '2) No'
            disp('opted out of reusing the existing centerline');
            delineation_flag = 1;
    end
else
    delineation_flag = 1;
end

%draw the centerline
if delineation_flag == 1
    disp('use image to map the centerline following flow...');
    beep;
    
    %trace the centerline
    disp('trace the centerline starting from the seaward side of the AOI, overlapping each end');
    figure(fig);
    % [~,~,~,xi,yi] = improfile; %old (SLOW) approach: be patient, a + cursor will appear, double right click to terminate
    pline = drawpolyline('Color','c','Linewidth',1); %a + cursor will appear, double right click to terminate
    xi = pline.Position(:,1); yi = pline.Position(:,2); clear pline;
    % plot(xi,yi,'xy'); hold on;
    
    %create a regularly interpolated version of the centerline
    prof_length = nansum(sqrt((xi(2:end)-xi(1:end-1)).^2 + (yi(2:end)-yi(1:end-1)).^2));
    interpCL_100m = interparc(round(prof_length/100),xi,yi,'spline'); %create a spooth line with spacing of ~100 m
    % plot(interpCL_100m(:,1),interpCL_100m(:,2),'sy'); hold on;
    prof_length = nansum(sqrt((interpCL_100m(2:end,1)-interpCL_100m(1:end-1,1)).^2 + (interpCL_100m(2:end,2)-interpCL_100m(1:end-1,2)).^2));
    interpCL = interparc(round(prof_length/spacer),xi,yi,'spline'); %create a spooth line with spacing defined by "spacer"
    C.X = interpCL(:,1); C.Y = interpCL(:,2);
    %plot
    plot(C.X,C.Y,'-y','linewidth',2); %plot the centerline
    drawnow;
    %clear out old variables
    clear xi yi interpCL* prof_length;
end

%convert the centerline coordinates to along-profile distance from the origin
CL(1) = 0;
for k = 2:length(C.X)
    CL(k) = CL(k-1)+sqrt((C.X(k)-C.X(k-1)).^2 + (C.Y(k)-C.Y(k-1)).^2);
    CL_ang(k-1) = atand((C.Y(k)-C.Y(k-1))./(C.X(k)-C.X(k-1)));
end
CL_ang(length(C.X)) = CL_ang(length(C.X)-1);


%% extend lines orthogonal from the centerline orientation at even increments
xy = [S.X, S.Y];

%calculate the distance from each centerline point to the nearest AOI
%vertex then use the maximum distance to define the length of the
%orthogonal vectors extending from the centerline
for j = 1:length(CL)
    dist = sqrt((C.X(j)-S.X).^2 + (C.Y(j)-S.Y).^2);
    min_dist(j) = min(dist,[],'all','omitnan');
    clear dist;
end
halfwidth = max(min_dist); %define the half-width vector
if halfwidth < 1000; halfwidth = 1000; end

%create transects
inc = transect_inc/spacer;
XF.X = []; XF.Y = []; ind = 1;
center_points = [1:inc:length(CL)];
%if the fjord polygon is really long, decrease transect spacing increment to avoid memory issues and crashes
if length(center_points) > 10 && halfwidth > 5000
    transect_inc = 1.5*transect_inc; 
    inc = transect_inc/spacer;
    clear center_points; center_points = [1:inc:length(CL)];
    transect_spacer = 2*spacer;
else
    transect_spacer = spacer;
end
%loop through all the centerline points & extend transects
for j = center_points
    disp([num2str(ceil(j/inc)),' of ',num2str(ceil(length(CL)/inc)),' transects']);

    %check that the transect's intersection with the centerline falls
    %within the fjord polygon
    [stat] = inpoly2([C.X(j),C.Y(j)],xy);
    if stat == 1
        %look for intersections with the AOI, incrementally increasing
        %length as needed
        k=1;
        while k
            l=0:transect_spacer:(1+(k*0.25))*halfwidth;

            %orthogonal +90 degrees from the centerline orientation
            perp = CL_ang(j)+90;
            del_y = l.*sind(perp); del_x = l.*cosd(perp);
            yp = C.Y(j)+del_y; xp = C.X(j)+del_x;
            %orthogonal -90 degrees from the centerline orientation
            perp_neg = CL_ang(j)-90;
            del_y = l.*sind(perp_neg); del_x = l.*cosd(perp_neg);
            yn = C.Y(j)+del_y; xn = C.X(j)+del_x;
            %combine the vectors
            xt = [fliplr(xn),xp(2:end)]; yt = [fliplr(yn),yp(2:end)];

            %find intersections of transect and fjord outline
            [xi,yi] = polyxpoly(S.X,S.Y,xt,yt);

            %check that there are intersections on both sides of the centerline
            dist_pts = sqrt((xi-C.X(j)).^2 + (yi-C.Y(j)).^2);
            dir_pts = sign(atan2d((yi-C.Y(j)),(xi-C.X(j))));
            vec_pts = dir_pts.*dist_pts;

            %too many intersections, filter!
            if length(xi) > 2 && ~isempty(find(vec_pts<0)) && ~isempty(find(vec_pts>0))
                %find the indices for closest intersecting points
                pt1_ind = find(vec_pts==max(vec_pts(vec_pts<0)));
                pt2_ind = find(vec_pts==min(vec_pts(vec_pts>0)));

                %isolate only the closest intersections
                xi_temp = xi; yi_temp = yi; clear xi yi;
                xi = xi_temp([pt1_ind,pt2_ind],1);
                yi = yi_temp([pt1_ind,pt2_ind],1);
                dist_pts = dist_pts([pt1_ind,pt2_ind]);
                dir_pts = dir_pts([pt1_ind,pt2_ind]);
                vec_pts = vec_pts([pt1_ind,pt2_ind]);
                clear pt*_ind *i_temp;

                %crop the transect so that the multiple intersections
                %don't lead to finding weird points inside the polygon
                %that make the outline incorrect
                dist_tran = sqrt((xt-C.X(j)).^2 + (yt-C.Y(j)).^2);
                dir_tran = sign(atan2d((yt-C.Y(j)),(xt-C.X(j))));
                vec_tran = dir_tran.*dist_tran;
                vec_inds = find(vec_tran>vec_pts(vec_pts<0) & vec_tran<vec_pts(vec_pts>0));
                xt = xt([min(vec_inds)-1,vec_inds,max(vec_inds)+1]); yt = yt([min(vec_inds)-1,vec_inds,max(vec_inds)+1]);

                clear d*_tran vec_tran;
            else %multiple intersections but only on one side so transect needs to be lengthened
                k = k + 1;
            end

            % crop the transect and save it if it intersects one point on
            % each side of the centerline
            if length(xi) == 2 && ~isempty(find(vec_pts<0)) && ~isempty(find(vec_pts>0))
                %find the transect points inside the AOI
                in = inpolygon(xt,yt,S.X,S.Y); inrefs = find(in==1);
                %add the intersection points to the appropriate ends of the
                %transect within the AOI
                edge_dist = sqrt((xt(in)-xi(1)).^2 + (yt(in)-yi(1)).^2);
                if find(edge_dist==min(edge_dist)) == 1
                    x_perp = [xi(1),xt(inrefs(1:end)),xi(2)]; y_perp = [yi(1),yt(inrefs(1:end)),yi(2)];
                elseif find(edge_dist==min(edge_dist)) == sum(in)
                    x_perp = [xi(1),xt(inrefs(end:-1:1)),xi(2)]; y_perp = [yi(1),yt(inrefs(end:-1:1)),yi(2)];
                else
                    error('something went wrong with transect cropping!');
                end
                clear in inrefs edge_dist;

                %add transects to a structure if they still intersect the centerline
                [xc,~] = polyxpoly(C.X,C.Y,x_perp,y_perp);
                if ~isempty(xc)
                    XF(ind).X = x_perp; XF(ind).Y = y_perp;
                    plot(x_perp,y_perp,'--y','linewidth',2); hold on; drawnow;
                    ind = ind+1;
                end
                clear xc x_perp y_perp;

                break;
            else %lengthen the transect
                k = k+1;
            end
        end
    else
        disp('transect falls outside the AOI');
    end
    clear stat d*_pts vec_pts;
    clear perp* del_* xp yp xn ynxt yt xi yi;
end

%% check for overlapping transects and temporarily modify the melange mask to block transect overlap

%identify overlaps
for j = 1:length(XF)
    for k = 1:length(XF)
        if k ~= j
            [xi,yi] = polyxpoly([XF(j).X(1),XF(j).X(end-1)],[XF(j).Y(1),XF(j).Y(end-1)],[XF(k).X(1),XF(k).X(end-1)],[XF(k).Y(1),XF(k).Y(end-1)]);
            if ~isempty(xi)
                transect_crosses(j,k) = k;
                [xic,yic] = polyxpoly(XF(j).X,XF(j).Y,C.X,C.Y);
                cross_dist(j,k) = sqrt((xi-xic).^2 + (yi-yic).^2);
                clear xic yix
            else
                transect_crosses(j,k) = NaN;
                cross_dist(j,k) = NaN;
            end
        else
            transect_crosses(j,k) = NaN;
            cross_dist(j,k) = NaN;
        end
        clear xi yi;
    end
end
%identify the coordinates of the crossing transects that are closest to the centerline
[mindist,minind] = min(cross_dist,[],'all');
[j,k] = ind2sub(size(cross_dist),minind);
[xi,yi] = polyxpoly([XF(j).X(1),XF(j).X(end-1)],[XF(j).Y(1),XF(j).Y(end-1)],[XF(k).X(1),XF(k).X(end-1)],[XF(k).Y(1),XF(k).Y(end-1)]);
split_ind = max([j,k]);

%create a temporary melange mask that is modified as needed to provide a
%barrier for otherwise overlapping transects
if sum(~isnan(transect_crosses),'all') > 0
    transect_overlap = 1; %flag that the transects overlap and the melange mask is modifed
    
    % %find the middle of neighboring intersections
    % for j = 1:length(XF)
    %     if ismember(j-1,transect_crosses(j,:)) || ismember(j+1,transect_crosses(j,:))
    %         %         neighbor_crosses(j) = sum([ismember(j-1,transect_crosses(j,:)),ismember(j+1,transect_crosses(j,:))]);
    %         neighbor_crosses(j) = 1;
    %     else
    %         neighbor_crosses(j) = 0;
    %     end
    % end
    % mid_crosses = round(nanmean(find(neighbor_crosses==1)));

    %split the intersecting lines to create the transect splitter
    %jth transect
    dxj = mode(diff(XF(j).X)); dyj = mode(diff(XF(j).Y));
    [xjc1,yjc1,ijc1] = polyxpoly([XF(j).X(1)-dxj,xi],[XF(j).Y(1)-dyj,yi],C.X,C.Y);
    [xjc2,yjc2,ijc2] = polyxpoly([XF(j).X(end-1)+2*dxj,xi],[XF(j).Y(end-1)+2*dyj,yi],C.X,C.Y);
    %correct vector extending to the polygon edge will not intersect the centerline
    if isempty(xjc1) 
        [xjp,yjp] = polyxpoly([XF(j).X(1)-dxj,xi],[XF(j).Y(1)-dyj,yi],S.X,S.Y); %intersect with the outline
        XF(j).X = [xi,XF(j).X(ijc2(1,1):end)]; XF(j).Y = [yi,XF(j).Y(ijc2(1,1):end)]; %now crop the transect
    else
        [xjp,yjp] = polyxpoly([xi,XF(j).X(end-1)+2*dxj],[yi,XF(j).Y(end-1)+2*dyj],S.X,S.Y); %intersect with the outline
        XF(j).X = [XF(j).X(1:ijc1(1,1)),xi]; XF(j).Y = [XF(j).Y(1:ijc1(1,1)),yi]; %now crop the transect
    end
    %kth transect
    dxk = mode(diff(XF(k).X)); dyk = mode(diff(XF(k).Y));
    [xkc1,ykc1,ikc1] = polyxpoly([XF(k).X(1)-dxk,xi],[XF(k).Y(1)-dyk,yi],C.X,C.Y);
    [xkc2,ykc2,ikc2] = polyxpoly([XF(k).X(end-1)+2*dxk,xi],[XF(k).Y(end-1)+2*dyk,yi],C.X,C.Y);
    %correct vector extending to the polygon edge will not intersect the centerline
    if isempty(xkc1) 
        [xkp,ykp] = polyxpoly([XF(k).X(1)-dxk,xi],[XF(k).Y(1)-dyk,yi],S.X,S.Y); %intersect with the outline
        XF(k).X = [xi:XF(k).X(ikc2(1,1):end)]; XF(k).Y = [yi:XF(k).Y(ikc2(1,1):end)]; %now crop the transect
    else
        [xkp,ykp] = polyxpoly([xi,XF(k).X(end-1)+2*dxk],[yi,XF(k).Y(end-1)+2*dyk],S.X,S.Y); %intersect with the outline
        XF(k).X = [XF(k).X(1:ikc1(1,1)),xi]; XF(k).Y = [XF(k).Y(1:ikc1(1,1)),yi]; %now crop the transect
    end
    %join the points then extend the lin to be sure it intersects the outline
    xp = [nanmean([xjp,xkp]),xi]; yp = [nanmean([yjp,ykp]),yi]; 
    dxp = diff(xp); dyp = diff(yp);
    xd = [xp-dxp*0.5, xp]; yd = [yp-dyp*0.5, yp]; 
    %find the intersections for the transect and the outline
    [xs,ys,is] = polyxpoly(xd,yd,S.X,S.Y);
    %create a temp melange mask & add the part of the line that has
    %intersections to the shape
    melpoly_x = S.X(1:is(1,2)); melpoly_y = S.Y(1:is(1,2));
    %only include the transect to the center-most intersection
    melpoly_x = [melpoly_x; xs(1); xi; xs(1); melmask.uncropped.x(is(1,2)+1:end)];
    melpoly_y = [melpoly_y; ys(1); yi; ys(1); melmask.uncropped.y(is(1,2)+1:end)];

    
    % %determine what side of the centerline the intersections occur on (assuming
    % %it is just one side!) then find the segment of the fjord polygon that the
    % %transect intersects on that side
    % j = mid_crosses;
    % for k = 1:length(XF)
    %     if k ~= j
    %         if ~isempty(polyxpoly([XF(j).X(1),XF(j).X(end-1)],[XF(j).Y(1),XF(j).Y(end-1)],[XF(k).X(1),XF(k).X(end-1)],[XF(k).Y(1),XF(k).Y(end-1)]))
    %             [xi(k),yi(k)] = polyxpoly([XF(j).X(1),XF(j).X(end-1)],[XF(j).Y(1),XF(j).Y(end-1)],[XF(k).X(1),XF(k).X(end-1)],[XF(k).Y(1),XF(k).Y(end-1)]);
    %         else
    %             xi(k) = NaN; yi(k) = NaN;
    %         end
    %     else
    %         xi(k) = NaN; yi(k) = NaN;
    %     end
    % end
    % %extend the transect so it intersects the melange outline
    % dx1 = mode(diff(XF(j).X)); dy1 = mode(diff(XF(j).Y));
    % if size(XF(j).X,1) == 1
    %     xd = [XF(j).X(1)-dx1,XF(j).X(1:end-1),XF(j).X(end-1)+dx1];
    %     yd = [XF(j).Y(1)-dy1,XF(j).Y(1:end-1),XF(j).Y(end-1)+dy1];
    % else
    %     xd = [XF(j).X(1)-dx1;XF(j).X(1:end-1);XF(j).X(end-1)+dx1];
    %     yd = [XF(j).Y(1)-dy1;XF(j).Y(1:end-1);XF(j).Y(end-1)+dy1];
    % end
    % %find the intersections for the transect and the outline
    % [xs,ys,is] = polyxpoly(xd,yd,melmask.uncropped.x,melmask.uncropped.y);
    % %calculate the distance between each transect intersection point with the
    % %middle of the intersecting transects and each melange edge
    % dists(1,:) = sqrt((xs(1)-xi).^2 + (ys(1)-yi).^2);
    % dists(2,:) = sqrt((xs(2)-xi).^2 + (ys(2)-yi).^2);
    % %find which melange edge is closer to the transect intersections
    % ind = find(dists == min(dists,[],'all')); [row,~] = ind2sub(size(dists),ind);
    % 
    % %create a temp melange mask & add the part of the line that has
    % %intersections to the shape
    % melpoly_x = melmask.uncropped.x(1:is(row,2)); melpoly_y = melmask.uncropped.y(1:is(row,2));
    % %only include the transect to the center-most intersection
    % melpoly_x = [melpoly_x; xs(row); xi(find(dists(row,:) == max(dists(row,:)))); xs(row); melmask.uncropped.x(is(row,2)+1:end)];
    % melpoly_y = [melpoly_y; ys(row); yi(find(dists(row,:) == max(dists(row,:)))); ys(row); melmask.uncropped.y(is(row,2)+1:end)];
    
else
    transect_overlap = 0;
    
    melpoly_x = S.X;
    melpoly_y = S.Y;
end
clear *_crosses xi yi xd yd xs ys is dists ind row dx* dy* x*c1 x*c2 y*c1 y*c2 xkp ykp *xp *yp;
plot(melpoly_x,melpoly_y,'-k','linewidth',2);

%now crop the transects
for j = 12:length(XF)
    %find the intersection with the centerline
    [xij,yij,ij] = polyxpoly([XF(j).X,XF(j).X(end-1)],[XF(j).Y,XF(j).Y(end-1)],AF.X,AF.Y);

    %find intersections with the melange outline
    dx2 = mode(diff(XF(j).X)); dy2 = mode(diff(XF(j).Y));
    if size(XF(j).X,1) == 1
        XF(j).X = [XF(j).X(1)-dx2,XF(j).X(1:end-1),XF(j).X(end-1)+dx2];
        XF(j).Y = [XF(j).Y(1)-dy2,XF(j).Y(1:end-1),XF(j).Y(end-1)+dy2];
    else
        XF(j).X = [XF(j).X(1)-dx2;XF(j).X(1:end-1);XF(j).X(end-1)+dx2];
        XF(j).Y = [XF(j).Y(1)-dy2;XF(j).Y(1:end-1);XF(j).Y(end-1)+dy2];
    end
    [xi,yi,ii] = polyxpoly(XF(j).X,XF(j).Y,melpoly_x,melpoly_y);

    %if there are >2 intersections because the fjord curved, filter out
    %the ones farthest from the centerline... the 'first' & 'last' are
    %somewhat arbitrary but prevent double identification of non-unique
    %intersections where a hinge in the polygons doubles-back
    if size(ii,1) > 2 %first transect
        di = ii(:,1)-ij(:,1);
        ref1a = find(di==max(di(di<0)),1,'last'); ref1b = find(di==min(di(di>0)),1,'first');
        XF(j).X = XF(j).X(ii(ref1a,1):ii(ref1b,1)); XF(j).Y = XF(j).Y(ii(ref1a,1):ii(ref1b,1));
        clear xi yi ii di ref1*;
    end

    %now check that the transects don't still intersect, cropping those
    %that come after the hinge as needed
    if j > split_ind
        for k = split_ind+1:j-1
            [xt,~] = polyxpoly([XF(j).X(1),XF(j).X(end-1)],[XF(j).Y(1),XF(j).Y(end-1)],[XF(k).X(1),XF(k).X(end-1)],[XF(k).Y(1),XF(k).Y(end-1)]);
            if ~isempty(xt)
                [xt,yt,it] = polyxpoly(XF(j).X,XF(j).Y,XF(k).X,XF(k).Y);
                %assume the smaller segment should be removed
                seg1_dist = sqrt((XF(j).X(1)-xt).^2 + (XF(j).Y(1)-yt).^2);
                seg2_dist = sqrt((XF(j).X(end-1)-xt).^2 + (XF(j).Y(end-1)-yt).^2);
                if seg1_dist > seg2_dist
                    XF(j).X = [XF(j).X(1:it(1,1)),xt]; XF(j).Y = [XF(j).Y(1:it(1,1)),yt];
                else
                    XF(j).X = [xt,XF(j).X(it(1,1):end)]; XF(j).Y = [yt,XF(j).Y(it(1,1):end)];
                end
            end
        end
    end


    clear xij yij ij dx2 dy2 xt yt seg*_dist;
    plot(XF(j).X,XF(j).Y,'-c','linewidth',2); hold on; drawnow;
end


%% export the polylines (as a csv or shapefile depending on selection)

%make a directory to house all shapefiles if it doesn't exist
if ~isfolder([site_dir,'shapefiles']); mkdir([site_dir,'shapefiles']); end

%assume that the last directory specified in the site_dir path is named
%after the study site & use it to define a site_abbrev variable used in the
%centerline and transect exported file names
slashes = strfind(site_dir,'/');
if slashes(end) == length(site_dir) %if site_dir has a trailing /
    site_abbrev = site_dir(slashes(find(slashes<length(site_dir),1,'last'))+1:length(site_dir)-1);
else %if site_dir does not have a trailing /
    site_abbrev = site_dir(slashes(end+1):length(site_dir));
end

if strcmpi(export_type,'shp')
    disp('exporting centerline & transects as shapefiles');
    %centerline
    s.Geometry = 'Polyline';
    s.BoundingBox = double([min(C.X) min(C.Y); max(C.X) max(C.Y)]);
    s.X = double(C.X); s.Y = double(C.Y);
    shapewrite(s,[site_dir,'shapefiles/',site_abbrev,'-centerline.shp']);
    s.CoordinateReferenceSystem = site_crs;
    wkt = wktstring(s.CoordinateReferenceSystem);
    writematrix(wkt,[site_dir,'shapefiles/',site_abbrev,'-centerline.prj'],'FileType','text', 'QuoteStrings', false);
    clear s;

    %transects
    for j = 1:length(XF)
        s(j).Geometry = 'Polyline';
        s(j).BoundingBox = double([min(S.X) min(S.Y); max(S.X) max(S.Y)]);
        s(j).X = double(XF(j).X); s(j).Y = double(XF(j).Y);
    end
    shapewrite(s,[site_dir,'shapefiles/',site_abbrev,'-transects_',num2str(transect_inc),'m.shp']);
    writematrix(wkt,[site_dir,'shapefiles/',site_abbrev,'-transects_',num2str(transect_inc),'m.prj'],'FileType','text', 'QuoteStrings', false);
    clear s;
else
    disp('exporting centerline & transects as CSVs');

    %centerline
    T=table(C.Y,C.X);
    column_names = ["Northing (m)","Easting (m)"];
    T.Properties.VariableNames = column_names;
    writetable(T,[site_dir,'shapefiles/',site_abbrev,'-centerline.csv']);
    clear T;

    %transects
    TY = []; TX = [];
    for j = 1:length(XF)
        TY = [TY; XF(j).Y'; NaN];
        TX = [TX; XF(j).X'; NaN];
    end
    T=table(TY,TX);
    column_names = ["Northing (m)","Easting (m)"];
    T.Properties.VariableNames = column_names;
    writetable(T,[site_dir,'shapefiles/',site_abbrev,'-transects_',num2str(transect_inc),'m.csv']);
    clear T TY TX;

end

%% create CSV of coordinates for the intersection of each transect with the centerline (for velocity extraction)
xy = [S.X, S.Y]; 
w = who;

%create centerline profile vector
if ismember('C',w) == 0 && ismember('AF',w) == 1
    C.X = AF.X; C.Y = AF.Y;
end

inc = transect_inc/spacer; ind = 1;
for j = 1:inc:length(C.X)
    %check that the transect's intersection with the centerline falls
    %within the fjord polygon, then add it if it does
    [stat] = inpoly2([C.X(j),C.Y(j)],xy);
    if stat == 1
        %add coordinates to the vector
        AX_X(ind,1) = C.X(j); AX_Y(ind,1) = C.Y(j);

        %advance
        ind = ind+1;
        clear stat;
    end
end
%throw an error if there are not the same number of centerline points as transects
if length(AX_X) ~= length(XF)
    error('Check code because centerline points & transects do not have the same number')
end

%write to file (uncomment the if/else statement if you'd prefer to default
%to outputs in the same format as the centerline & transect profiles)
% if strcmpi(export_type,'shp')
    disp('exporting centerline points at transect intersections as shapefiles');
    %centerline
    s.Geometry = 'PolyLine';
    s.BoundingBox = double([min(AX_X) min(AX_X); max(AX_Y) max(AX_Y)]);
    s.X = double(AX_X); s.Y = double(AX_Y);
    shapewrite(s,[site_dir,'shapefiles/',site_abbrev,'-centerline_',num2str(transect_inc),'m-interval.shp']);
    s.CoordinateReferenceSystem = site_crs;
    wkt = wktstring(s.CoordinateReferenceSystem);
    writematrix(wkt,[site_dir,'shapefiles/',site_abbrev,'-centerline_',num2str(transect_inc),'m-interval.prj'],'FileType','text', 'QuoteStrings', false);
    clear s;

% else
    disp('exporting centerline points at transect intersections as CSVs');

    %centerline
    AX_X(isnan(AX_X)) = []; AX_Y(isnan(AX_Y)) = []; 
    if size(AX_X,2) == 1
        T=table(AX_Y,AX_X);
    else
        T=table(AX_Y',AX_X');
    end
    column_names = ["Northing (m)","Easting (m)"];
    T.Properties.VariableNames = column_names;
    writetable(T,[site_dir,'shapefiles/',site_abbrev,'-centerline_',num2str(transect_inc),'m-interval.csv']);
    clear T;

% end

end