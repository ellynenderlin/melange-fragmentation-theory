% export melange mask shapefiles
clearvars; close all;

%root directory for all glacier folders
root_dir = '/Users/ellynenderlin/Research/NSF_Greenland-Calving/fragmentation-curves/';
EPSG_file = '/Users/ellynenderlin/Research/miscellaneous/EPSG3413.prj';

%identify the site folders
cd_to_root = ['cd ',root_dir]; eval(cd_to_root);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.')
        sitenames = [sitenames; sites(i).name];
    end
end

%loop through the folders & extract info
disp('Exporting melange bounding boxes...');
for i = 1:length(sitenames)
    cd_to_site = ['cd ',sitenames(i,:),'/']; eval(cd_to_site);

    %load the matfile containing the melange masks
    load_file = ['load ',sitenames(i,:),'-melange-masks.mat']; eval(load_file);
    S.Geometry = 'Polygon';
    S.BoundingBox = [max(melmask.uncropped.x) min(melmask.uncropped.y);
        min(melmask.uncropped.x) max(melmask.uncropped.y)];
    S.X = melmask.uncropped.x; S.Y = melmask.uncropped.y;
    S.Name = sitenames(i,:);
    S.Projection = 'EPSG:3413';
    shapewrite(S,[sitenames(i,:),'-melange-outline.shp']); 
    copyfile(EPSG_file,[root_dir,sitenames(i,:),'/',sitenames(i,:),'-melange-outline.prj']);
    
    %clear and continue loop
    clear S;
    cd ..
end
    
    