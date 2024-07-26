function [melmask,velSA,velSAmad] = analyze_melange_velocity_coherence(site_abbrev,root_dir,outline_dir,outline_file,melmask)

%%% Examine melange velocity time series to determine when and where there is coherent melange.
% addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
%     '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/');

%specify directories & files
% site_abbrev = 'ASG'; %3-letter abbrevation
% root_dir = '/Volumes/Jokulhaup_5T/Greenland-melange/';
site_dir = [root_dir,site_abbrev,'/']; %overarching directory for sites
vel_dir = [site_dir,'velocities/']; %site-specific directory holding ITS_LIVE time series csvs for each point
% outline_dir = site_dir; outline_file = [site_abbrev,'-melange-masks.mat'];
cd(site_dir);

%specify years of interest (from Sentinel 2 launch to near-present)
vel_years = [2016:1:2022]; %ITS_LIVE coverage only through 2022 as of June 1st, 2024

%only complete the velocity-coherence analysis if it has not been done
if ~isfield(melmask.dated,'Vext_x')
    %load all the data: number csv files first so that farthest from the
    %glacier is 01 & numbers increase approaching the glacier
    vels = dir([vel_dir,'*.csv']); warning off;
    for j = 1:size(vels,1)
        %read in the data for the point
        T = readtable([vel_dir,vels(j).name]);
        V(j).lat = T.lat; V(j).lon = T.lon;
        [V(j).x,V(j).y] = wgs2ps(V(j).lon,V(j).lat);
        V(j).date = T.mid_date; V(j).dt = T.date_dt_days_;
        V(j).v = T.v_m_yr_; V(j).v_err = T.v_error_m_yr_;
        clear T;
    end
    clear vels;
    
    % %plot the data for a first look
    % cmap = colormap(parula(3));
    % figure;
    % for j = 1:length(V)
    %     scatter(V(j).date,V(j).v,10,cmap(j,:)); hold on;
    % end
    
    %find the extent of the year when there is a continuous data record to use
    %as an estimate of when the melange covers each point
    for j = 1:length(V)
        %create a vector of the date range for each velocity estimate
        datenum_mid = datenum(V(j).date);
        datenum_range = [datenum_mid - (V(j).dt/2), datenum_mid + (V(j).dt/2)];
        datestr_mid = datestr(V(j).date,'yyyymmdd');
        datestr_min = datestr(datenum_range(:,1),'yyyymmdd'); datestr_max = datestr(datenum_range(:,2),'yyyymmdd');
        for k = 1:size(datenum_mid,1)
            decidate_mid(k,:) = convert_to_decimaldate(datestr_mid(k,:));
            decidate_min(k,:) = convert_to_decimaldate(datestr_min(k,:));
            decidate_max(k,:) = convert_to_decimaldate(datestr_max(k,:));
        end
        [sorted_dates,daterefs] = sort(decidate_mid);
        sorted_decidate_min = decidate_min(daterefs); sorted_decidate_max = decidate_max(daterefs);
        
        %     %plot the data to get a first look
        %     figure;
        %     errorbar(decidate_mid,V(j).v,V(j).v_err,V(j).v_err,abs(decidate_mid-decidate_min),abs(decidate_mid-decidate_max),'sb');
        
        %find periods of continuous data coverage for each year
        for k = 1:length(vel_years)
            year_refs(k) = numel(find(sorted_dates>=vel_years(k) & sorted_dates<vel_years(k)+1));
            first_ref = find(sorted_decidate_min>=vel_years(k),1,'first');
            start_ref(k) = first_ref; start_date(k) = sorted_decidate_min(start_ref(k));
            p=1;
            while p
                next_ref = find(sorted_dates>sorted_dates(first_ref),1,'first');
                if sorted_decidate_min(next_ref) <= sorted_decidate_max(first_ref) %date ranges overlap
                    first_ref = next_ref; clear next_ref;
                    p=p+1;
                elseif sorted_dates(next_ref) - sorted_dates(first_ref) < 30/365 %speeds are less than a month apart
                    first_ref = next_ref; clear next_ref;
                    p=p+1;
                else %decent gap in speed from melange attributed to lack of melange coherence (aka melange break-up)
                    end_ref(k) = first_ref; end_date(k) = sorted_decidate_max(first_ref);
                    end_datestr(k,:) = datestr_max(find(decidate_max == end_date(k),1,'last'),:);
                    break
                end
            end
            clear first_ref;
            
            %replace persistent melange (or on-ice!) observations with 1 (1=full year)
            if end_date(k) > vel_years(k)+1
                end_date(k) =  vel_years(k)+1;
                end_datestr(k,:) = [num2str(vel_years(k)),'1231'];
            end
        end
        
        %add to velocity data structure
        %     V(j).mel_loss_yrs = vel_years(year_refs>=50);
        %     V(j).mel_loss_decidate = end_date(year_refs>=50);
        %     V(j).mel_loss_datestr = end_datestr(year_refs>=50,:);
        V(j).mel_loss_yrs = vel_years;
        V(j).mel_loss_decidate = end_date;
        V(j).mel_loss_datestr = end_datestr;
        
        clear datenum_* datestr_* *decidate* sorted_* daterefs* year_refs start_* end_*;
    end
    
    %save the compiled velocity data
    save([vel_dir,site_abbrev,'-melange-velocities.mat'],'V','-v7.3');
    
    %plot the melange break-up dates
    % figure1 = figure; %plot by year
    figure2 = figure; %plot by point
    yr_cmap = colormap(cool(length(vel_years)));
    for k = 1:length(vel_years)
        for j = 1:length(V)
            if ~isempty(find(V(j).mel_loss_yrs-vel_years(k)==0))
                mel_loss_decidate(j,k) = V(j).mel_loss_decidate(find(V(j).mel_loss_yrs-vel_years(k)==0))-vel_years(k);
            else
                mel_loss_decidate(j,k) = NaN;
            end
        end
        mel_loss_decidate(mel_loss_decidate==0) = NaN;
        %     figure(figure1);
        %     plot(mel_loss_decidate(:,k)-vel_years(k),[1:1:length(V)],'-x','color',yr_cmap(k,:)); hold on;
    end
    %make a boxplot of the timing for each point
    figure(figure2);
    boxchart((mel_loss_decidate)'); hold on;
    disp('Check that the break-up dates (fraction of year) in the boxplot seem reasonable before moving on');
    
    % convert typical melange extent for each season to an area estimate using the velocity time series & melange shapefile
    close all;
    if exist('V') ~= 1; load([vel_dir,site_abbrev,'-melange-velocities.mat']); end
    
    %load the melange shapefile that was used to download DEMs
    load([outline_dir,outline_file]);
    S.X = melmask.uncropped.x; S.Y = melmask.uncropped.y;
    
    %check which velocity points are inside the melange polygon
    for j = 1:length(V); point_flag(j) = inpolygon(V(j).x(1),V(j).y(1),S.X,S.Y); end
    
    %plot the shape and the velocity points to double-check
    figure;
    plot([S.X;S.X(1)],[S.Y;S.Y(1)],'-k','linewidth',2); axis xy equal; hold on;
    for j = 1:length(V); if point_flag(j) == 1; plot(V(j).x(1),V(j).y(1),'xr'); hold on; end; end
    drawnow;
    
    %draw long perpendicular bisectors from each velocity point
    for j = 1:length(V)
        line_length = 10000; %line length in m (10km)
        if j == 1
            line_angle(j,:) = 90+atan2d((V(j+1).x(1)-V(j).x(1)),(V(j+1).y(1)-V(j).y(1))); %flip x & y because I want the angle 90 degrees from the centerline
        elseif j > 1 && j < length(V)
            line_angle(j,:) = 90+atan2d(mean([(V(j+1).x(1)-V(j).x(1)),(V(j).x(1)-V(j-1).x(1))]),mean([(V(j+1).y(1)-V(j).y(1)),(V(j).y(1)-V(j-1).y(1))])); %flip x & y because I want the angle 90 degrees from the centerline
        else
            line_angle(j,:) = 90+atan2d((V(j).x(1)-V(j-1).x(1)),(V(j).y(1)-V(j-1).y(1))); %flip x & y because I want the angle 90 degrees from the centerline
        end
        line_y(j,:) = [V(j).y(1)+line_length*cosd(line_angle(j,:)), V(j).y(1)-line_length*cosd(line_angle(j,:))];
        line_x(j,:) = [V(j).x(1)+line_length*sind(line_angle(j,:)), V(j).x(1)-line_length*sind(line_angle(j,:))];
    end
    
    %use the median of the end image date used to create velocities over the
    %specified time period to identify when melange breaks up
    modays_norm = [31 28 31 30 31 30 31 31 30 31 30 31];
    cumdays_norm = cumsum(modays_norm); cumdays_norm = [0 cumdays_norm];
    mel_loss_decidate_med = nanmedian(mel_loss_decidate,2);
    for k = 1:12
        lineref = find(mel_loss_decidate_med>=(cumdays_norm(k)/max(cumdays_norm)),1,'first');
        if isempty(lineref)
            mel_loss_ref(k) = 0;
        else
            mel_loss_ref(k) = lineref;
        end
        clear lineref;
    end
    
    %for each DEM date, use the seasonal down-glacier melange centerline extent
    %and the manually-delineated terminus to crop the melange shapefile and
    %solve for the melange surface area
    clear S; date_cmap = colormap(jet(length(melmask.dated)));
    for j = 1:length(melmask.dated)
        %load the melange mask that has been cropped to the terminus
        S.X = melmask.dated(j).x; S.Y = melmask.dated(j).y;
        berg_mo(j,:) = str2num(melmask.dated(j).datestring(5:6));
        
        %find the corresponding seaward extent for the month
        sea_ref = mel_loss_ref(str2num(melmask.dated(j).datestring(5:6))); %seaward extent reference based on velocities
        disp(['j = ',num2str(j),' & vel ref = ',num2str(sea_ref)]);
        if sea_ref ~= 1
            %check that the melange bisector point is not behind the terminus
            in = inpolygon(mean(line_x(sea_ref,:)),mean(line_y(sea_ref,:)),melmask.dated(j).x,melmask.dated(j).y);
            if sum(in) == 0
                melmask.dated(j).Vext_x = melmask.dated(j).x; melmask.dated(j).Vext_y =melmask.dated(j).y;
                disp(['No coherent melange for ',melmask.dated(j).datestring,' according to velocities']);
            else
                %find where the melange bisector intersects the fjord mask
                out_intercept = []; out_interceptx = []; out_intercepty = [];
                for i = 1:length(S.X)-1
                    [xi,yi] = polyxpoly(S.X(i:i+1),S.Y(i:i+1),line_x(sea_ref,:),line_y(sea_ref,:)); %find the intersections of the terminus trace with the melange outline
                    if ~isempty(xi)
                        out_intercept = [out_intercept i]; out_interceptx = [out_interceptx xi]; out_intercepty = [out_intercepty yi];
                    end
                    clear xi yi;
                end
                
                %create a new melange polygon cropped on both ends
                melmask.dated(j).Vext_x = [out_interceptx(find(out_intercept==min(out_intercept))); S.X(min(out_intercept)+1:max(out_intercept)-1); out_interceptx(find(out_intercept==max(out_intercept)));out_interceptx(find(out_intercept==min(out_intercept)))];
                melmask.dated(j).Vext_y = [out_intercepty(find(out_intercept==min(out_intercept))); S.Y(min(out_intercept)+1:max(out_intercept)-1); out_intercepty(find(out_intercept==max(out_intercept)));out_intercepty(find(out_intercept==min(out_intercept)))];
                
                clear out_intercept*;
            end
        else
            melmask.dated(j).Vext_x = melmask.dated(j).x;
            melmask.dated(j).Vext_y = melmask.dated(j).y;
        end
        plot(melmask.dated(j).Vext_x,melmask.dated(j).Vext_y,'-','color',date_cmap(j,:),'linewidth',2); axis xy equal; hold on;
        
        clear S sea_ref;
    end
    disp('Melange masks re-saved with velocity-based extent estimates.')
    save([outline_dir,outline_file],'melmask','-v7.3');
else
    %identify month of each DEM 
    for j = 1:length(melmask.dated)
        berg_mo(j,:) = str2num(melmask.dated(j).datestring(5:6));
    end
end

%solve for seasonal surface areas
DJF = find(berg_mo==12 | berg_mo == 1 | berg_mo == 2);
MAM = find(berg_mo==3 | berg_mo == 4 | berg_mo == 5);
JJA = find(berg_mo==6 | berg_mo == 7 | berg_mo == 8);
SON = find(berg_mo==9 | berg_mo == 10 | berg_mo == 11);
%calculate median and MAD surface areas for each season
disp('Solving for seasonal velocity coherence-based area estimates...');
for k = 1:length(DJF); melpoly = polyshape(melmask.dated(DJF(k)).x,melmask.dated(DJF(k)).y); SA(k) = area(melpoly); clear melpoly; end
velSA(1) = nanmedian(SA); velSAmad(1) = mad(SA,1); clear SA;
for k = 1:length(MAM); melpoly = polyshape(melmask.dated(MAM(k)).x,melmask.dated(MAM(k)).y); SA(k) = area(melpoly); clear melpoly; end
velSA(2) = nanmedian(SA); velSAmad(2) = mad(SA,1); clear SA;
for k = 1:length(JJA); melpoly = polyshape(melmask.dated(JJA(k)).x,melmask.dated(JJA(k)).y); SA(k) = area(melpoly); clear melpoly; end
velSA(3) = nanmedian(SA); velSAmad(3) = mad(SA,1); clear SA;
for k = 1:length(SON); melpoly = polyshape(melmask.dated(SON(k)).x,melmask.dated(SON(k)).y); SA(k) = area(melpoly); clear melpoly; end
velSA(4) = nanmedian(SA); velSAmad(4) = mad(SA,1); clear SA;
disp('Typical seasonal melange areas estimated from velocities');

end