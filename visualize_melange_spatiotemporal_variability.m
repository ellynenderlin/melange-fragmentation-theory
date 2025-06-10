%%% Create plots & maps of GrIS melange data

clearvars; close all; warning off;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/inpoly2/');
addpath('/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-fragmentation-code/');

%specify paths
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-melt/';
years = 2011:1:2023;

%identify the site folders
cd(root_dir);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.')
        sitenames = [sitenames; sites(i).name];
    end
end

%loop through the folders & extract info
disp('Average elevation of each transect plotted over time & centerline velocity');
MP = struct;
for j = 1:length(sitenames)
    disp(sitenames(j,:));
    cd([root_dir,sitenames(j,:),'/']);
    MP(j).name = sitenames(j,:);

    %open the CSV of time-stamped elevation transects
    T = readtable([sitenames(j,:),'_transects_elevations.csv']);
    T_headers = T.Properties.VariableNames;
    datestart = strfind(string(T_headers(3)),'20');
    for p = 3:size(T,2)
        MP(j).date(1,p-2) = string(T_headers{p}(datestart:datestart+7));
    end

    %find the NaNs in the coordinate pairs to identify each transect
    coords = table2array(T(:,1:2));
    nan_inds = find(isnan(coords(:,1))==1);
    
    %extract the data
    tran_ind = 1; k = 1; zprofs = [];
    while k < size(T,1)
        if k == 1
            %extract full coordinates for transect
            MP(j).XF(tran_ind).X = coords(1:nan_inds(k)-1,2);
            MP(j).XF(tran_ind).Y = coords(1:nan_inds(k)-1,1);

            %extract the centroid coordinates
            MP(j).XF(tran_ind).centroid = nanmean(table2array(T(1:nan_inds(k)-1,1:2)));

            %calculate the median elevation for each date
            MP(j).XF(tran_ind).Z = nanmean(table2array(T(1:nan_inds(k)-1,3:end)));

            %add to a temp matrix for plotting
            zprofs(tran_ind,:) = MP(j).XF(tran_ind).Z;

            tran_ind = tran_ind+1;
        elseif ismember(k,nan_inds)
            if ~ismember(k+1,nan_inds) %check for back-to-back NaNs
                %extract full coordinates for transect
                MP(j).XF(tran_ind).X = coords(k+1:nan_inds(find(nan_inds==k)+1)-1,2);
                MP(j).XF(tran_ind).Y = coords(k+1:nan_inds(find(nan_inds==k)+1)-1,1);

                %extract the centroid coordinates
                MP(j).XF(tran_ind).centroid = nanmean(table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,1:2)));

                %calculate the median elevation for each date
                MP(j).XF(tran_ind).Z = nanmean(table2array(T(k+1:nan_inds(find(nan_inds==k)+1)-1,3:end)));

                %add to a temp matrix for plotting
                zprofs(tran_ind,:) = MP(j).XF(tran_ind).Z;

                tran_ind = tran_ind+1;
            end
        end
        k = k+1;
    end
    clear T T_headers datestart coords nan_inds;

    %plot time-series of the width-averaged elevation profiles
    figure; set(gcf,'position',[50 50 500 1000]);
    sub_yr = subplot(2,1,1); yr_cmap = cmocean('matter',length(years)+2); yr_cmap = yr_cmap(2:end,:);
    sub_mo = subplot(2,1,2); mo_cmap = cmocean('phase',12);
    mean_prof = nanmean(zprofs,2); end_ind = find(~isnan(mean_prof)==1,1,'first');
    for p = 1:size(zprofs,2)
        yr = str2num(MP(j).date{p}(1:4)); mo = str2num(MP(j).date{p}(5:6));
        subplot(sub_yr);
        plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',zprofs(end_ind:end,p),'color',yr_cmap(yr-min(years)+1,:),'linewidth',2); hold on; 
        %add transect width plotted on the secondary y-axis
        
        subplot(sub_mo);
        plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',zprofs(end_ind:end,p),'color',mo_cmap(mo,:),'linewidth',2); hold on;

    end
    subplot(sub_yr);
    title(sitenames(j,:))
    plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',mean_prof(end_ind:end),'color','k','linewidth',1.5); hold on; 
    grid on; drawnow;
    subplot(sub_mo); 
    plot([0:2000:2000*(size(zprofs,1)-1-(end_ind-1))]',mean_prof(end_ind:end),'color','k','linewidth',1.5); hold on; 
    grid on; drawnow;

    clear zprofs *_cmap mean_prof end_ind;
end








