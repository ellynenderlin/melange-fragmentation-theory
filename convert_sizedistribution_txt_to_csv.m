%%% Convert size distribution textfiles to csvs 
clearvars; close all;

%navigate to overarching directory containing site sub-directories
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange-melt/';
cd(root_dir);

%specify the column headers for the txt & csv files
column_names = ["Count", "SurfaceArea_mean", "SurfaceArea_range"];
column_units = ["icebergs", "m^2", "m^2"];

%go into each site sub-directory & re-write tab-delimited textfiles as csvs with column names
sites = dir([root_dir,'*']);
for j = 1:size(sites,1)
    if exist(sites(j).name) == 7 && length(sites(j).name) == 3 %3-letter abbreviation folder
        cd(sites(j).name); disp(sites(j).name);
        
        %set up dummy vectors to compile data for the site
        berg_nos = []; berg_A = []; berg_dA = []; berg_dates = []; berg_dateformat = [];
        
        %loop through date-specific textfiles
        berg_dists = dir('*-iceberg-distribution.txt');
        if isempty(berg_dists)
            disp('no iceberg distribution textfiles in directory');
        else
            for k = 1:size(berg_dists,1)
                %load the date-specific textfile
                txt = dlmread(berg_dists(k).name);
                txt(txt(:,1)==0,1) = NaN;
                
                %add data to the site-compiled table
                if isempty(berg_A)
                    berg_A = [txt(:,2)];
                    berg_dA = [txt(:,3)];
                end
                berg_nos = [berg_nos, txt(:,1)];
%                 berg_dates = [berg_dates; repmat(convert_to_decimaldate(berg_dists(k).name(5:12)),size(txt(:,1)))];
                berg_dates = [berg_dates, string(berg_dists(k).name(5:12))];
                berg_dateformat = [berg_dateformat, "YYYYMMDD"];
                
                %re-export the date-specific size distributions as a csv
                T=table(txt(:,1),txt(:,2),txt(:,3));
                T.Properties.VariableNames = column_names; T.Properties.VariableUnits = column_units;
                writetable(T,[berg_dists(k).name(1:end-4),'.csv']);
                
                %clear dated info
                clear txt T;
            end
            clear berg_dists;
            
            %export the data compiled for the site as a single table
%             berg_nos(berg_nos == 0) = NaN;
            T=array2table([berg_A,berg_dA,berg_nos]);
            site_column_names = ["SurfaceArea_mean", "SurfaceArea_range",berg_dates];
            site_column_units = ["m^2", "m^2", berg_dateformat];
            T.Properties.VariableNames = site_column_names; T.Properties.VariableUnits = site_column_units;
            writetable(T,[sites(j).name,'-melange-distributions.csv']);
            clear T;
        end
        
        %move to the next site
        cd ..
    end
end
