%%% This code was writen by Karina Zikan to statistically compare elevations
%%% from multiple ICESat-2 products to reference elevations calculated by
%%% DTM_reference_elevations_calculation.m
%%%
%%% SPECIFIED INPUTS:
%%%     icesat2 = path to ICESat-2 file (output from batch_icesat2_terrain_parameter_comparison.m)
%%%     ref_elevations = path to reference elevations file (output from DTM_reference_elevations_calculation.m)
%%%     abbrev = site abriviation for file name
%%%     snowcover = set to snow-on ('snowon') or snow-off ('snowoff') conditions
%%% OUTPUTS:
%%%     
%%%
%%% Last updated: May 2023 by Karina Zikan

%% Inputs
clearvars;
addpath(['./functions'])
addpath(['/Users/karinazikan/Documents/cmocean'])

%Folder path 
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/RCEW/';
%site abbreviation for file names
abbrev = 'RCEW';
%Set snowcover to 'snowonn' or 'snowoff'
snowcover = 'snowoff';
%Turn slope correction off or on
slope_correction = 0; % 0 = off, 1 = on

%% Set file paths and colormaps
%File paths
icesat2_atl08 = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL08-params'];
ref_elevations_atl08 = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL08-ref-elevations'];

icesat2_atl06 = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr'];
ref_elevations_atl06 = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-params'];

icesat2_atl06_class = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-atl08class'];
ref_elevations_atl06_class = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-atl08class-params'];

%set colors
colors{1} = cmocean('-dense',6);
colors{2} = cmocean('-algae',5);
colors{3} = cmocean('ice',5);

%% Load data
%load the reference elevation data
E_08 = readtable(ref_elevations_atl08);
E_06 = readtable(ref_elevations_atl06);
E_06_class = readtable(ref_elevations_atl06_class);

%load the ICESat-2 data
I_08 = readtable(icesat2_atl08); %read in files
I_06 = readtable(icesat2_atl06);
I_06_class = readtable(icesat2_atl06_class);

footwidth = 11; % approx. width of icesat2 shot footprint in meters

%% Filter snow-on or snow-off for ATL08
% % % % % bright = I_08.Brightness_Flag;
% % % % % dates_08 = datetime(I_08.date);

% season (1=winter(Jan,Feb,Mar),2=spring(Apr,May,Jun),3=summer(Jul,Aug,Sep),4=fall(Oct,Nov,Dec))
if snowcover == 'snowoff'
    ib = find(I_08.season == 3); 
    disp('Snow off')
elseif snowcover == 'snowonn'
    ib = find(I_08.season == 1);
    disp('Snow on')
else
    error('snowcover must be set to snowonn or snowoff')
end
I_08 = I_08(ib,:);
E_08 = E_08(ib,:);

% Filter snow-on or snow-off for ATl06
dates_06 = datetime(I_06.time.Year,I_06.time.Month,I_06.time.Day);
if snowcover == 'snowoff'
    ib = find(I_06.time.Month >=6 & I_06.time.Month <= 9);
    disp('Snow off')
elseif snowcover == 'snowonn'
    ib = find(I_06.time.Month <=2 | I_06.time.Month >=12);
    disp('Snow on')
else
    error('snowcover must be set to snowon or snowoff')
end
I_06 = I_06(ib,:);
E_06 = E_06(ib,:);

% Filter snow-on or snow-off for ATl06 atl08class
dates_06_class = datetime(I_06_class.time.Year,I_06_class.time.Month,I_06_class.time.Day);
if snowcover == 'snowoff'
    ib = find(I_06_class.time.Month >=6 & I_06_class.time.Month <= 9);
    disp('Snow off')
elseif snowcover == 'snowonn'
    ib = find(I_06_class.time.Month <=2 | I_06_class.time.Month >=12);
    disp('Snow on')
else
    error('snowcover must be set to snowonn or snowoff')
end
I_06_class = I_06_class(ib,:);
E_06_class = E_06_class(ib,:);


%% Loop through each product
N_Products = 3; %number of products to analyze
for i = 1:N_Products
    if i == 1
        T = I_08;
        E = E_08;
        acronym = 'ATL08';

        % ICESat-2 data
        zmod = T.Elevation_bestfit(:); % save the fitted 'model' elevations (icesat-2 elevations)
        zmod = T.Elevation(:); % save the 'model' elevations (icesat-2 elevations)
        zstd = T.std; %save the standard deviation of the icesat-2 elevation estimates
        easts = T.Easting(:); % pull out the easting values
        norths = T.Northing(:); % pull out the northings
        slope = T.slope(:);
        aspect = T.aspect(:);
        %canopy = T.Canopy(:);
    elseif i == 2
         T = I_06;
        E = E_06;
        acronym = 'ATL06sr';    

        % ICESat-2 data
        zmod = T.h_mean(:); % save the median 'model' elevations (icesat-2 elevations
        zstd = T.h_sigma; %save the standard deviation of the icesat-2 elevation estimates
        easts = T.Easting(:); % pull out the easting values
        norths = T.Northing(:); % pull out the northings
        slope = E.slope_mean(:);
        aspect = E.aspect_mean(:);
%        canopy = T.Canopy(:);
    else
        T = I_06_class;
        E = E_06_class;
        acronym = 'ATL06sr atl08 classifications';
        % ICESat-2 data
        zmod = T.h_mean(:); % save the median 'model' elevations (icesat-2 elevations
        zstd = T.h_sigma; %save the standard deviation of the icesat-2 elevation estimates
        easts = T.Easting(:); % pull out the easting values
        norths = T.Northing(:); % pull out the northings
        slope = E.slope_mean(:);
        aspect = E.aspect_mean(:);
     %   canopy = T.Canopy(:);
    end
    % Ref elev data
    elevation_report(:,1) = E.elevation_report_nw_mean;  %non-weighted mean ref elevation
    elevation_report(:,2) = E.elevation_report_mean;    %weighted mean ref elevation
    elevation_report(:,3) = E.elevation_report_fitted;    %weighted & fitted ref elevation
    elevation_report_std = E.elevation_report_std;  %std of elevations within footprint

    %% Stats
    %calculate the elevation residuals
    for j = 1:3
        differences = zmod - elevation_report(:,j); %calculate the icesat2 elevations and the calculated reference elevations
        differences(differences > 80) = NaN; differences(differences < -80) = NaN; %remove extreme outliers
        Dmean{i}(:,j) = nanmean(differences); % calculate mean of diferences
        Dstd{i}(:,j) = std(differences,'omitnan'); % calculate std of diferences
        Dmed{i}(:,j) = median(differences,'omitnan');
        Dmad{i}(:,j) = median(abs(differences-Dmean{i}(:,j)),'omitnan');
        zrmse{i}(:,j) = sqrt(nansum((differences).^2)./length(differences)); %calculate rmse of  differeces

        %        Dks_test(i,:) = kstest(differences); %kolmagorov smirnof test

        %         % Removing residuals below -13
        %         ix = find(differences < -13);
        %         zmod(ix) = NaN;
        %         differences(ix) = NaN;

        % Vertically corregister
        if slope_correction == 0
            % Vertical corregistration
            Residuals(:,j) = differences-Dmed{i}(:,j); 
        elseif slope_correction == 1
            Residuals = differences-Dmed{i}(:,j);
            %calculate quadratic slope correction
            x= slope; y = Residuals;
            ind = isnan(x) | isnan(y); %index nans
            x(ind) = []; y(ind) = []; %remove nans
            p = polyfit(x,y,2); % fit quadratic
            % Vertical corregistration
            Residuals(:,j) = differences-Dmed{i}(:,j)-polyval(p,slope);
        else
            error('slope_correction must be set to 0 (no slope correction) or 1 (slope correction applied)')
        end
        %   Residuals(:,1) is comparison to non weighted mean elevations
        %   Residuals(:,2) is comparison to weighted mean elevations
        %   Residuals(:,3) is comparison to weighted & fitted elevations

        Rmean{i}(:,j) = nanmean(Residuals(:,j)); % calculate mean of Residuals
        Rstd{i}(:,j) = std(Residuals(:,j),'omitnan'); % calculate std of Residuals
        Rmed{i}(:,j) = median(Residuals(:,j),'omitnan');
        Rmad{i}(:,j) = median(abs(Residuals(:,j)-Rmean{i}(:,j)),'omitnan');
        Rrmse{i}(:,j) = sqrt(nansum((Residuals(:,j)).^2)./length(Residuals(:,j))); %calculate rmse of  Residuals
        

    end
    Residuals(:,4) = i;
    ResidualsAll{i} = Residuals;
    
    %% Plots
    % Reference historgrams
    fig2 = figure(2);
    subplot(N_Products,1,i); set(gcf,'position',[50 50 800 500]); clear h;
    binwidth = 0.2;
    h(1) = histogram(Residuals(:,1),'Normalization','pdf'); h(1).BinWidth = binwidth; h(1).FaceAlpha = 1; h(1).FaceColor = colors{i}(1,:);  h(1).EdgeColor = 'k'; hold on;
    h(2) = histogram(Residuals(:,2),'Normalization','pdf'); h(2).BinWidth = binwidth; h(2).FaceAlpha = 0.75; h(2).FaceColor = colors{i}(2,:); h(2).EdgeColor = 'k';
    h(3) = histogram(Residuals(:,3),'Normalization','pdf'); h(3).BinWidth = binwidth; h(3).FaceAlpha = 0.75;  h(3).FaceColor = colors{i}(3,:); h(3).EdgeColor = 'w';
    plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
    set(gca,'fontsize',16,'xlim',[-4 2]);
    legend(h,['DEM: Non-weighted mean, std = ' num2str(Dstd{i}(:,1))],['DEM: Weighted mean, std = ' num2str(Dstd{i}(:,2))],['DEM: Weighted and fitted, std = ' num2str(Dstd{i}(:,3))],'Location','northwest');
    xlabel('Vertical offset (m)'); ylabel('Probability density'); title(acronym);

    clear elevation_report Residuals
end
%% Plots outside loop
% Products historgrams
% fig4 = figure(4); clf; hold on
% for i = 1:N_Products
% %     if i == 1
% %         h(i) = histogram(ResidualsAll{i}(:,3),'Normalization','pdf');  h(i).FaceAlpha = 1; h(i).BinWidth = binwidth; h(i).FaceColor = colors{i}(3,:);  h(i).EdgeColor = 'k';
% %         fplot(@(x) mynormpdf(x,nanmean(ResidualsAll{i}(:,3)), std(ResidualsAll{i}(:,3),'omitnan')),[-10 8], 'Linewidth', 2,'Color',colors{i}(2,:));
% %     else
%         h(i) = histogram(ResidualsAll{i}(:,1),'Normalization','pdf');  h(i).FaceAlpha = .5; h(i).BinWidth = binwidth; h(i).FaceColor = colors{i}(3,:);  h(i).EdgeColor = 'k';
%         fplot(@(x) mynormpdf(x,nanmean(ResidualsAll{i}(:,1)), std(ResidualsAll{i}(:,1),'omitnan')),[-10 8], 'Linewidth', 2,'Color',colors{i}(2,:));
% %     end
% end
% %plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
% set(gca,'fontsize',16,'xlim',[-4 2]);
% set(gcf,'position',[50 50 800 400]);
% legend('ATL08','ATL08 pdf','ATL06sr','ATL06sr pdf','ATL06sr ATL08 clasifications','ATL06sr ATL08 clasifications pdf');
% xlabel('Vertical offset (m)'); ylabel('Probability density');
% txt = {['N-08 = ' num2str(length(ResidualsAll{1}(:,1))-sum(isnan(ResidualsAll{1}(:,1))))],['N-06 = ' num2str(length(ResidualsAll{2}(:,1))-sum(isnan(ResidualsAll{2}(:,1))))],['N-06_{class} = ' num2str(length(ResidualsAll{3}(:,1))-sum(isnan(ResidualsAll{3}(:,1))))]};
% text(-8,.2,txt);

%Non-parametirc pdfs w/ histograms
fig5 = figure(5); clf; hold on
for i = 1:N_Products
%     if i == 1
%         h(i) = histogram(ResidualsAll{i}(:,3),'Normalization','pdf');  h(i).FaceAlpha = 1; h(i).BinWidth = binwidth; h(i).FaceColor = colors{i}(3,:);  h(i).EdgeColor = 'k';
%         pd = fitdist(ResidualsAll{i}(:,3),'kernel','Kernel','normal'); 
%         fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(2,:));
%     else
        h(i) = histogram(ResidualsAll{i}(:,1),'Normalization','pdf');  h(i).FaceAlpha = .5; h(i).BinWidth = binwidth; h(i).FaceColor = colors{i}(3,:);  h(i).EdgeColor = 'k';
        pd = fitdist(ResidualsAll{i}(:,1),'kernel','Kernel','normal'); 
        fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(2,:));
%     end
end
%plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
set(gca,'fontsize',16,'xlim',[-4 2]);
set(gcf,'position',[50 50 800 400]);
legend('ATL08','ATL08 pdf','ATL06sr','ATL06sr pdf','ATL06sr ATL08 class','ATL06sr ATL08 class pdf');
xlabel('Vertical offset (m)'); ylabel('Probability density');
txt = {['N-08 = ' num2str(length(ResidualsAll{1}(:,1))-sum(isnan(ResidualsAll{1}(:,1))))],['N-06 = ' num2str(length(ResidualsAll{2}(:,1))-sum(isnan(ResidualsAll{2}(:,1))))],['N-06_{class} = ' num2str(length(ResidualsAll{3}(:,1))-sum(isnan(ResidualsAll{3}(:,1))))]};
text(-3,.2,txt);

fig6 = figure(6); clf; hold on
for i = 1:N_Products
%     if i == 1
%         pd = fitdist(ResidualsAll{i}(:,3),'kernel','Kernel','normal'); 
%         fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
%     else
        pd = fitdist(ResidualsAll{i}(:,1),'kernel','Kernel','normal'); 
        fplot(@(x) pdf(pd,x),[-4 4], 'Linewidth', 3,'Color',colors{i}(3,:));
%     end
end
set(gca,'fontsize',16);
set(gcf,'position',[50 50 800 400]);
legend('ATL08','ATL06','Classified ATL06','Location','northwest');
xlabel('Vertical offset (m)'); ylabel('Probability density');
hold off

% % Normal pdfs for each tested product
% fig7 = figure(7); clf; hold on
% for i = 1:N_Products
% %     if i == 1
% %         fplot(@(x) mynormpdf(x,Dmean{i}(:,1), Dstd{i}(:,1)),[-10 8], 'Linewidth', 3,'Color',colors{i}(3,:));
% %     else
%         fplot(@(x) mynormpdf(x,Dmean{i}(:,1), Dstd{i}(:,1)),[-10 8], 'Linewidth', 3,'Color',colors{i}(3,:));
% %     end
% end
% set(gca,'fontsize',16);
% set(gcf,'position',[50 50 800 400]);
% legend('ATL08','ATL06sr','ATL06sr ATL08 clasifications','Location','northwest');
% xlabel('Vertical offset (m)'); ylabel('Probability density');
% hold off

%% create boxcharts for each terrain parameter
ResidualsTable.residuals = [ResidualsAll{1}(:,3); ResidualsAll{2}(:,1); ResidualsAll{3}(:,1)];
ResidualsTable.product = [ResidualsAll{1}(:,4); ResidualsAll{2}(:,4); ResidualsAll{3}(:,4)];
ResidualsTable.elevations = [E_08.elevation_report_mean; E_06.elevation_report_mean; E_06_class.elevation_report_mean];
ResidualsTable.aspect = [I_08.aspect; E_06.aspect_mean; E_06_class.aspect_mean];
ResidualsTable.slope = [I_08.slope; E_06.slope_mean; E_06_class.slope_mean];

% ResidualsTable.residuals = [ResidualsAll{2}(:,1); ResidualsAll{3}(:,1)];
% ResidualsTable.product = [ResidualsAll{2}(:,4); ResidualsAll{3}(:,4)];
% ResidualsTable.elevations = [E_06.elevation_report_mean; E_06_class.elevation_report_mean];
% ResidualsTable.aspect = [I_06.aspect; I_06_class.aspect];
% ResidualsTable.slope = [I_06.slope; I_06_class.slope];

Nbin = 8; %set humber of bins for box plots
fig8 = figure(8); clf
%ELEVATION
subplot(3,1,1);
h = histogram(E_08.elevation_report_mean(~isnan(ResidualsAll{1}(:,1))),Nbin);
elev_binwidth = h.BinWidth; elev_binedges = h.BinEdges;
clear h;
xlabel('Elevation (m a.s.l.)','fontsize',16); %ylabel('Elevation residuals (m)','fontsize',16);
%close(gcf);
%ASPECT
subplot(3,1,2);
h = histogram(I_08.aspect(~isnan(ResidualsAll{1}(:,1))),Nbin);
aspect_binwidth = h.BinWidth; aspect_binedges = h.BinEdges;
clear h;
xlabel('Aspect (degrees)','fontsize',16); ylabel('Observations','fontsize',16);
% close(gcf);
%SLOPE
subplot(3,1,3);
h = histogram(I_08.slope(~isnan(ResidualsAll{1}(:,1))),Nbin);
slope_binwidth = h.BinWidth; slope_binedges = h.BinEdges;
clear h;
xlabel('Slope (degrees)','fontsize',16);
% close(gcf);

whiskerline = '-'; outliermarker = 'o';

% % boxplot figure
% fig9 = figure(9); clf
% %ELEVATION
% subplot(3,1,1);
% hold on
% bins = {num2str(elev_binedges(2)) num2str(elev_binedges(3)) num2str(elev_binedges(4)) num2str(elev_binedges(5)) num2str(elev_binedges(6))};
% groupElev = discretize(ResidualsTable.elevations,elev_binedges,'categorical',bins);
% %meanWeight = groupsummary(ResidualsTable.residuals,groupElev,'mean');
% boxchart(groupElev,ResidualsTable.residuals,'GroupByColor',ResidualsTable.product,'MarkerStyle','none')
% colororder([colors{1}(3,:); colors{2}(3,:); colors{3}(3,:)]); 
% ylim([-5,3])
% % plot(meanWeight,'-o','Color','k')
% set(gca,'fontsize',16,'box','on'); drawnow;
% legend('ATL08','ATL06sr','ATL06sr ATL08 clasifications','Location','northwest');
% xlabel('Elevation (m a.s.l.)','fontsize',16); %ylabel('Elevation residuals (m)','fontsize',16);
% %text(1,max(ylims)-0.05*range(ylims),'a)','fontsize',16);
% %ASPECT
% subplot(3,1,2);
% hold on
% bins = {num2str(aspect_binedges(2)) num2str(aspect_binedges(3)) num2str(aspect_binedges(4)) num2str(aspect_binedges(5)) num2str(aspect_binedges(6))};
% groupAspect = discretize(ResidualsTable.aspect,aspect_binedges,'categorical',bins);
% %meanWeight = groupsummary(ResidualsTable.residuals,groupAspect,'mean');
% boxchart(groupAspect,ResidualsTable.residuals,'GroupByColor',ResidualsTable.product,'MarkerStyle','none')
% ylim([-5,3])
% % plot(meanWeight,'-o','Color','k')
% set(gca,'fontsize',16,'box','on'); drawnow;
% xlabel('Aspect (degrees)','fontsize',16); ylabel('Elevation residuals (m)','fontsize',16);
% %text(1,max(ylims)-0.05*range(ylims),'a)','fontsize',16);
% %SLOPE
% subplot(3,1,3);
% hold on
% bins = {num2str(slope_binedges(2)) num2str(slope_binedges(3)) num2str(slope_binedges(4)) num2str(slope_binedges(5)) num2str(slope_binedges(6))};
% groupSlope = discretize(ResidualsTable.aspect,slope_binedges,'categorical',bins);
% %meanWeight = groupsummary(ResidualsTable.residuals,groupAspect,'mean');
% boxchart(groupSlope,ResidualsTable.residuals,'GroupByColor',ResidualsTable.product,'MarkerStyle','none')
% ylim([-5,3])
% % plot(meanWeight,'-o','Color','k')
% set(gca,'fontsize',16,'box','on'); drawnow;
% xlabel('Slope (degrees)','fontsize',16); %ylabel('Elevation residuals (m)','fontsize',16);
% %text(1,max(ylims)-0.05*range(ylims),'a)','fontsize',16);

% boxplot figure seperate products
fig10 = figure(10); clf
%ELEVATION
bins = {num2str(elev_binedges(2))};
for i= 3:length(elev_binedges)
    bins = [bins; {num2str(elev_binedges(i))}];
end
% ATL08
subplot(3,3,1);
hold on
groupElev = discretize(E_08.elevation_report_mean,elev_binedges,'categorical',bins);
boxchart(groupElev,ResidualsAll{1}(:,3),'BoxFaceColor',colors{1}(3,:),'MarkerStyle','none')
ylim([-2,2])
set(gca,'fontsize',16,'box','on'); drawnow;
title('ATL08');
xlabel('Elevation (m a.s.l.)','fontsize',16); ylabel('Elevation residuals (m)','fontsize',16);
% ATL06
subplot(3,3,2);
hold on
groupElev = discretize(E_06.elevation_report_mean,elev_binedges,'categorical',bins);
boxchart(groupElev,ResidualsAll{2}(:,1),'BoxFaceColor',colors{2}(3,:),'MarkerStyle','none')
ylim([-2,2])
set(gca,'fontsize',16,'box','on'); drawnow;
title('ATL06');
xlabel('Elevation (m a.s.l.)','fontsize',16); 
% ATL06-class
subplot(3,3,3);
hold on
groupElev = discretize(E_06_class.elevation_report_mean,elev_binedges,'categorical',bins);
boxchart(groupElev,ResidualsAll{3}(:,1),'BoxFaceColor',colors{3}(3,:),'MarkerStyle','none')
ylim([-2,2])
set(gca,'fontsize',16,'box','on'); drawnow;
title('ATL06 classified');
xlabel('Elevation (m a.s.l.)','fontsize',16); 
%ASPECT
bins = {num2str(aspect_binedges(2))};
for i= 3:length(aspect_binedges)
    bins = [bins; {num2str(aspect_binedges(i))}];
end
% ATL08
subplot(3,3,4);
hold on
groupAspect = discretize(I_08.aspect,aspect_binedges,'categorical',bins);
boxchart(groupAspect,ResidualsAll{1}(:,3),'BoxFaceColor',colors{1}(3,:),'MarkerStyle','none')
ylim([-2,2])
set(gca,'fontsize',16,'box','on'); drawnow;
xlabel('Aspect (degrees)','fontsize',16); ylabel('Elevation residuals (m)','fontsize',16);
% ATL06
subplot(3,3,5);
hold on
groupAspect = discretize(E_06.aspect_mean,aspect_binedges,'categorical',bins);
boxchart(groupAspect,ResidualsAll{2}(:,1),'BoxFaceColor',colors{2}(3,:),'MarkerStyle','none')
ylim([-2,2])
set(gca,'fontsize',16,'box','on'); drawnow;
xlabel('Aspect (degrees)','fontsize',16); 
% ATL06-class
subplot(3,3,6);
hold on
groupAspect = discretize(E_06_class.aspect_mean,aspect_binedges,'categorical',bins);
boxchart(groupAspect,ResidualsAll{3}(:,1),'BoxFaceColor',colors{3}(3,:),'MarkerStyle','none')
ylim([-2,2])
set(gca,'fontsize',16,'box','on'); drawnow;
xlabel('Aspect (degrees)','fontsize',16); 
%SLOPE
bins = {num2str(slope_binedges(2))};
for i= 3:length(slope_binedges)
    bins = [bins; {num2str(slope_binedges(i))}];
end
% ATL08
subplot(3,3,7);
hold on
groupSlope = discretize(I_08.slope,slope_binedges,'categorical',bins);
boxchart(groupSlope,ResidualsAll{1}(:,3),'BoxFaceColor',colors{1}(3,:),'MarkerStyle','none')
ylim([-6,6])
set(gca,'fontsize',16,'box','on'); drawnow;
xlabel('Slope (degrees)','fontsize',16); ylabel('Elevation residuals (m)','fontsize',16);
% ATL06
subplot(3,3,8);
hold on
groupSlope = discretize(E_06.slope_mean,slope_binedges,'categorical',bins);
boxchart(groupSlope,ResidualsAll{2}(:,1),'BoxFaceColor',colors{2}(3,:),'MarkerStyle','none')
ylim([-3,3])
set(gca,'fontsize',16,'box','on'); drawnow;
xlabel('Slope (degrees)','fontsize',16); 
% ATL06-class
subplot(3,3,9);
hold on
groupSlope = discretize(E_06_class.slope_mean,slope_binedges,'categorical',bins);
boxchart(groupSlope,ResidualsAll{3}(:,1),'BoxFaceColor',colors{3}(3,:),'MarkerStyle','none')
ylim([-3,3])
set(gca,'fontsize',16,'box','on'); drawnow;
xlabel('Slope (degrees)','fontsize',16); 

%%
% boxplot figure seperate products
fig1 = figure(1); clf
for i = 1:3
    subplot(4,1,i);
    boxchart(ResidualsAll{i},'Notch','on','MarkerStyle','none')
end
subplot(4,1,4)
clear grp
grp(1:length(ResidualsAll{1}(:,1))) = 1; 
grp([length(grp)]:[length(grp)+length(ResidualsAll{2}(:,1))]) = 2;
grp([length(grp)]:[length(grp)+length(ResidualsAll{3}(:,1))]) = 3;
boxplot([ResidualsAll{1}(:,1);ResidualsAll{2}(:,1);ResidualsAll{3}(:,1)],grp,'Notch','on')


%% Save figs
% saveas(fig1,['/Users/karinazikan/Documents/figures/' abbrev '_atl08_tracks'],'png')
% saveas(fig2,['/Users/karinazikan/Documents/figures/' abbrev 'DEMmethods_hist'],'png')
% saveas(fig3,['/Users/karinazikan/Documents/figures/' abbrev 'DEMmethods_pdf'],'png')
% saveas(fig4,['/Users/karinazikan/Documents/figures/' abbrev 'prod_hists_wpdf'],'png')
% saveas(fig5,['/Users/karinazikan/Documents/figures/' abbrev 'terrain_boxplots'],'png')
% saveas(fig6,['/Users/karinazikan/Documents/figures/' abbrev 'prod_pdf'],'png')
