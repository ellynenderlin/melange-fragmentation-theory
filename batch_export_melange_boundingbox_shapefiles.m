function batch_export_melange_boundingbox_shapefiles(root_dir,EPSG_file)
% export melange mask shapefiles
close all;

%identify the site folders
cd(root_dir);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.')
        sitenames = [sitenames; sites(i).name];
    end
end

%loop through the folders & extract info
disp('Exporting melange bounding boxes...');
for i = 1:length(sitenames)
    cd([sitenames(i,:),'/']);

    %load the matfile containing the melange masks
    melangefile = dir('*-melange-masks.mat');
    load(melangefile(1).name);
    S.Geometry = 'Polygon';
    S.BoundingBox = [max(melmask.uncropped.x) min(melmask.uncropped.y);
        min(melmask.uncropped.x) max(melmask.uncropped.y)];
    S.X = melmask.uncropped.x; S.Y = melmask.uncropped.y;
    S.Name = sitenames(i,:);
    S.Projection = 'EPSG:3413';
    
    %save in shapefiles directory
    if exist('shapefiles') == 0
        mkdir('shapefiles')
    end
    shapewrite(S,[root_dir,sitenames(i,:),'/shapefiles/',melangefile(1).name(1:end-10),'-outline.shp']); 
    copyfile(EPSG_file,[root_dir,sitenames(i,:),'/shapefiles/',melangefile(1).name(1:end-10),'-outline.prj']);
    
    %clear and continue loop
    clear S;
    cd ..
end
    
end