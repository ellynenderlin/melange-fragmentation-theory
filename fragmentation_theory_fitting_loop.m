function [alpha,c1,c2,c3,c4,data_lims,errors] = fragmentation_theory_fitting_loop(fnames,iceberg_dates)

%%This script fits the fragmentation theory function to the iceberg
%size distribution data n(v). The function is of the form:
% n(v) = c1*v^{-a}*exp(-v/c2) + c3*dv*exp(-v/c4)
%and we are looking for parameters c = [c1 c2 a c3 c4]'
%given data n,v,dv.
%The script uses a brute force approach to searching, which could be
%improved with steepest descent search instead
%
%Author: Michal A. Kopera
%         Department of Mathematics
%         Boise State University
%Date:    05/12/2020


%initialize
clf; %tic

%choose error norm to use:
% 2 - L2 norm (default)
% Inf - maximum norm
norm_type = Inf;

%choose bin cutoff size to account for submarine melt-induced deviation
%from a power law for the smallest icebergs
berg_cutoff = 1;

%prepare a figure
% fitplot = figure; 
for p = 1:size(fnames,1); yrs(p) = str2num(fnames(p,4:7)); end
% set(fitplot,'position',[50 50 500 ceil(length(unique(yrs))/2)*300]);
%pull indices for each year
% time_cmap = colormap(parula(length(fnames)));
season_cmap = 0.75*[0.000	0.000	0.000; 0.251	0.000	0.294; 0.463	0.165	0.514;
0.600	0.439	0.671; 0.761	0.647	0.812; 0.906	0.831	0.910;
0.851	0.941	0.827; 0.651	0.859	0.627; 0.353	0.682	0.380; 
0.106	0.471	0.216; 0.000	0.267	0.106; 0.50 0.50 0.50]; close gcf;

% read data
pl = []; date_names = []; ref = 0;
nfiles = size(fnames,1);

for p = 1:nfiles

    load_data = ['load ',fnames(p,:)]; eval(load_data);
    disp(['Looping through iceberg distribution #',num2str(p),' of ',num2str(nfiles)]);
    v1 = double(m.melange.Asurfs(m.melange.bergs~=0)'); dv1 = double(m.melange.binwidth(m.melange.bergs~=0)'); 
    n1 = double(m.melange.bergs(m.melange.bergs~=0)')./dv1; clear m;
    
    %establish plot variables
    if p == find(yrs == str2num(fnames(p,4:7)),1,'first'); ref = 0; date_names = []; figure; end
    ref = ref+1; date_names = [date_names; fnames(p,4:11)];
    
    if sum((n1.*dv1).*v1) >= 175e6
        
        %remove the very small iceberg fragments of BIG icebergs since those
        %are simply elevation peaks in the largest icebergs
        v1 = v1(n1>0.01); n1 = n1(n1>0.01);
        
%         %remove the first two point(influenced by submarine melting)
%         %we may need a more automatic way of choosing how many points to
%         %neglect
%         v = v1(berg_cutoff:end); n = n1(berg_cutoff:end); 
%         v = v(~isnan(n)); n = n(~isnan(n)); 
        
        v = v1(~isnan(n1)); n = n1(~isnan(n1)); % NO REMOVAL OF FIRST POINTS
        
        %begin the search loop - compute error for every combination of i_cut, j_cut
        err_array = 1e12*ones(length(v)); %pre-allocate error array
        c_array = zeros(length(v),length(v),5);
        a_array = zeros(length(v),length(v));
        
        for i=4:length(v)-3
            for j=4:length(v)-3
                %create a function which computes error for parameter a
                %given the cutoff values i_cut and j_cut
                %fun = @(a) fit(v,n,i_cut(i),j_cut(j),a,norm_type);
                
                %find for which a the error is minimal
                %a = fminbnd(fun,1.5, 2.0);
                a = 2;
                %find parameters and error for the fit with new a
                [error, c] = fit(v,n,i,j,a,norm_type);
                
                err_array(i,j) = error;
                c_array(i,j,:) = c;
                a_array(i,j) = a;
            end
        end
        
        %find global minimum in the i_cut, j_cut landscape
        %basically looking for which cutoff values we get minimal error
        [min_val,idx]=min(err_array(:));
        
        [i,j]=ind2sub(size(err_array),idx);
        c = c_array(i,j,:);
        a = a_array(i,j);
        
        %rename variables
        alpha(p)=c(3);
        c1(p) = c(1); c2(p) = c(2); c3(p) = c(4); c4(p) = c(5);
        data_lims(p,:) = [i j]; errors(p) = min_val;
        
        %display results
        fprintf(fnames(p,4:11)); fprintf("\n");
        fprintf("[c1, c2, c3, c4] = [%e %e %e %e]\n",c(1),c(2),c(4),c(5));
        fprintf("alpha = %f\n",c(3));
        fprintf("==========================================\n")
        %     toc
        
        %plot result
        figure(p); set(gca,'FontSize',16);
        loglog(v1,n1,'+','MarkerSize',6,'LineWidth',1,'Color',season_cmap(str2num(fnames(p,8:9)),:)); hold on;
        loglog(v1,model(c,v1),'LineWidth',2,'Color',season_cmap(str2num(fnames(p,8:9)),:)); hold on;
        xlabel('surface area [m^2]'); ylabel('count');
%         if p == find(yrs == str2num(fnames(p,4:7)),1,'last')
%             legend(pl,date_names); 
%             if isinf(norm_type)
%                 saveas(gcf,[fnames(p,1:2),'-iceberg-size-distributions_',fnames(p,4:7),'-',num2str(berg_cutoff),'mcutoff-INFnorm.png'],'png'); 
%             else
%                 saveas(gcf,[fnames(p,1:2),'-iceberg-size-distributions_',fnames(p,4:7),'-',num2str(berg_cutoff),'mcutoff-L2norm.png'],'png'); 
%             end
%         end
        grid on; drawnow;
        % hold off
    else
        fprintf("==========================================\n")
        fprintf(fnames(p,4:11)); fprintf("\n");
        fprintf("DEM covers too small of an area to get reliable data\n")
        fprintf("==========================================\n")
        
        %plot just the data
        figure(p); set(gca,'FontSize',16);
        loglog(v1,n1,'+','MarkerSize',6,'LineWidth',1,'Color',season_cmap(str2num(fnames(p,8:9)),:)); hold on;
%         if p == find(yrs == str2num(fnames(p,4:7)),1,'last')
%             legend(pl,date_names); 
%             if isinf(norm_type)
%                 saveas(gcf,[fnames(p,1:2),'-iceberg-size-distributions_',fnames(p,4:7),'-',num2str(berg_cutoff),'mcutoff-INFnorm.png'],'png'); 
%             else
%                 saveas(gcf,[fnames(p,1:2),'-iceberg-size-distributions_',fnames(p,4:7),'-',num2str(berg_cutoff),'mcutoff-L2norm.png'],'png'); 
%             end
%         end
        grid on; drawnow;
        
        %fill-in NaNs as outputs
        alpha(p)=NaN;
        c1(p) = NaN; c2(p) = NaN; c3(p) = NaN; c4(p) = NaN;
        data_lims(p,:) = [NaN NaN]; errors(p) = NaN;
    end
    %clear variables
    clearvars -except nfiles berg_cutoff norm_type fitplot time_cmap season_cmap fnames pl ref iceberg_dates date_names a_option yrs alpha c1 c2 c3 c4 data_lims errors;
end

end

%%Functions

function n = model(c,v) % Eq. 1 both sides

c1 = c(1); c2 = c(2); a = c(3); c3 = c(4); c4 = c(5);
n = c1 * v.^(-a) .* exp(-v/c2) + c3 * exp(-v/c4); 

end

function [error,c] = fit(v,n,i_cut,j_cut,a,norm_type)
%the actual fitting function


%exponential fit - left side
vl = v(1:i_cut);
nl = n(1:i_cut);

m = length(vl);
A = [ones(m,1),vl];
lhs = log(nl) - log(vl.^(-a)); % linearize by taking hte log
b = (A'*A)\(A'*lhs);
cl(1) = exp(b(1));
cl(2) = -1/b(2);
cl(3) = a;

%exponential fit - right side
vr = v(j_cut:end);
nr = n(j_cut:end);

m = length(vr);
A = [ones(m,1),vr];
rhs = log(nr); % linearize by taking the log
b = (A'*A)\(A'*rhs);

cr(1) = exp(b(1));
cr(2) = -1/b(2);

%now use the linear approach as an initial guess to fit the entire function
c0 = [cl, cr]; %Initial guess
%c0(find(c0<0)) = 0;
c0 = abs(c0); %making sure initial guess is at least positive in all coeficients

opts = optimset('Display','off');
c = lsqcurvefit(@model, c0, v, n,[],[],opts);
error = norm((model(c,v)-n)./n,norm_type);

end

end