function [Z,data_mask,gap_mask] = sl_correction(Z, WV_sigma)
%% This function performs the sea level corrections to DEMs used for 
% automated extraction of iceberg size distributions.

%Author: Jukes Liu, Boise State University
%Date:    12/29/2020
disp('Starting sea level correction')
%create a data mask to remove the edges around the DEM with no data
data_mask = zeros(size(Z.z.ortho));
for i = 1:size(Z.z.ortho)
    if ~isempty(find(~isnan(Z.z.ortho(i,:))==1,1,'first'))
        data_range(1) = find(~isnan(Z.z.ortho(i,:))==1,1,'first');
        data_range(2) = find(~isnan(Z.z.ortho(i,:))==1,1,'last');
        data_mask(i,data_range(1):data_range(2)) = 1;
    end
    clear data_range;
end
%find big holes in the melange portion of the DEM
gap_mask = zeros(size(Z.z.ortho));
gap_mask(isnan(Z.z.ortho)) = 1; gap_mask(data_mask==0) = 0;
gap_mask = gap_mask.*Z.fjord.DEM_mask;
%remove holes <3 pixels from the mask since they should be
%relatively fine to interpolate across
stats = regionprops(logical(gap_mask),'PixelIdxList','Area');
for i = 1:length(stats)
    if stats(i).Area<9 %3x3 pixels
        gap_mask(stats(i).PixelIdxList) = 0;
    end
end

%first guess at sea level from minimum elevations: used to
%estimate the maximum iceberg size, which is used as the search
%window size for another pass at the sea level estimate
fjord_elevs = Z.fjord.DEM_mask.*Z.z.ortho; fjord_elevs(Z.fjord.DEM_mask==0) = NaN;
wdim = 501; whalf = (wdim-1)/2;
for j = 1+whalf:whalf:size(Z.z.ortho,1)-(whalf+1)
    y_sub((j-1)/whalf) = Z.y(j);
    for k = [1:whalf:size(Z.z.ortho,2),size(Z.z.ortho,2)]
        xmin = k-whalf; xmax = k+whalf;
        if xmin<0; xmin = 1; end
        if xmax>size(Z.z.ortho,2); xmax = size(Z.z.ortho,2); end
        z_sub = fjord_elevs(j-whalf:j+whalf,xmin:xmax);
        k_index = round((k+whalf-1)/whalf);
        x_sub(k_index) = Z.x(k);
        if nansum(nansum(~isnan(z_sub)))>0
            zmin((j-1)/whalf,k_index) = nanmean(min(z_sub));
        else
            zmin((j-1)/whalf,k_index) = 0;
        end
        clear z_sub;
    end
end
zmin(zmin==0) = nanmean(zmin(zmin>0));
%smooth sea level bias & extrapolate to edges of the domain
[x_slgrid,y_slgrid] = meshgrid(x_sub,y_sub);
[x_grid,y_grid] = meshgrid(Z.x,Z.y);
zmin_interp1 = interp2(x_slgrid,y_slgrid,zmin,x_grid,y_grid);
zmin_fjord1 = Z.fjord.DEM_mask.*zmin_interp1; zmin_fjord1(zmin_fjord1 == 0) = NaN;
%apply a threshold elevation to pick-out bergs from sea ice & bergy bits
z_adjusted = Z.z.ortho - zmin_fjord1; z_adjusted(z_adjusted<0) = NaN;
fjord_logical = Z.fjord.DEM_mask.*z_adjusted;
fjord_logical(fjord_logical< 3*WV_sigma) = 0; fjord_logical(fjord_logical>0) = 1; fjord_logical(isnan(fjord_logical)) = 0;
%use regionprops to extract iceberg thickness (median, max, min) assuming hydrostatic equilibrium, length & width, and area
bergL = regionprops(logical(fjord_logical),'MajorAxisLength');
berg_length = [];
for k = 1:length(bergL)
    berg_length = [berg_length bergL(k).MajorAxisLength];
end
max_dim = round(max(berg_length));
if max_dim > 1000; max_dim = 1001; end
if mod(max_dim,2)==0; max_dim = max_dim+1; end
whalf = ((max_dim + 100)-1)/2;

%make sure minimum elevations (which should be sea level) are
%extracted from regions at least the size of the largest
%individual iceberg (half-width = whalf)
disp('Calculating sea level offset');
k=0; increment = 50;
zfjord_filled = fjord_elevs;
zfjord_filled(isnan(zfjord_filled)) = median(zmin_fjord1(~isnan(zmin_fjord1))); %fill gaps with the median of the minimum elevations
zfjord_filled = Z.fjord.DEM_mask.*zfjord_filled; zfjord_filled(~Z.fjord.DEM_mask) = NaN;
zfjord_filled(zfjord_filled==0) = NaN;
for i = 1:increment:size(Z.z.ortho,1)
    ymin = i-whalf; ymax = i+whalf;
    y_index = round((i-1)/increment)+1;
    if ymin<1; ymin = 1; end
    if ymax>=size(Z.z.ortho,1); ymax = size(Z.z.ortho,1); end
    for j = 1:increment:size(Z.z.ortho,2)
        xmin = j-whalf; xmax = j+whalf;
        if xmin<1; xmin = 1; end
        if xmax>=size(Z.z.ortho,2); xmax = size(Z.z.ortho,2); end
        
        %downsample the melange elevations & calculate stats
        z_sub = zfjord_filled(ymin:ymax,xmin:xmax);
        x_index = round((j-1)/increment)+1;
        if nansum(nansum(~isnan(z_sub))) == 0;
            zmin_ds(y_index,x_index) = NaN;
        else
            z_median = nanmedian(min(z_sub));
            z_mad = mad(min(z_sub),1);
            z_sub(z_sub < (z_median-2*(1.4826*z_mad))) = NaN; z_sub(z_sub > (z_median+2*(1.4826*z_mad))) = NaN;
            zmin_ds(y_index,x_index) = nanmean(min(z_sub));
        end
        
        %identify which areas in the melange elevation map have elevation data
        z_sub = fjord_elevs(ymin:ymax,xmin:xmax);
        if nansum(nansum(~isnan(z_sub))) == 0
            fjordz_ds(y_index,x_index) = NaN;
        else
            fjordz_ds(y_index,x_index) = 1;
        end
        clear xmin xmax z_sub;
    end
    clear ymin ymax;
end
%remove outlier elevations for sea level
zmin_ds(zmin_ds>(nanmedian(zmin_ds(~isnan(fjordz_ds)))+6*1.4826*mad(zmin_ds(~isnan(fjordz_ds)),1)))=NaN;
zmin_ds(zmin_ds<(nanmedian(zmin_ds(~isnan(fjordz_ds)))-6*1.4826*mad(zmin_ds(~isnan(fjordz_ds)),1)))=NaN;
%remove sea level bias
[x_slgrid,y_slgrid] = meshgrid(Z.x(1:increment:size(Z.z.ortho,2)),Z.y(1:increment:size(Z.z.ortho,1)));
[x_grid,y_grid] = meshgrid(Z.x,Z.y);
zmin_interp = interp2(x_slgrid,y_slgrid,zmin_ds,x_grid,y_grid);
zmin_fjord = Z.fjord.DEM_mask.*zmin_interp; zmin_fjord(zmin_fjord==0) = NaN;
Z.z.sealevel = single(zmin_fjord);
clear zmin* zmin_* zfjord fjord_elevs fjordz_ds;
Z.z.adjusted = nansum(cat(3,Z.z.ortho,-Z.z.sealevel),3); Z.z.adjusted(Z.fjord.DEM_mask==0) = Z.z.ortho(Z.fjord.DEM_mask==0);
disp('Sea level correction finished')
end