%%% This code was writen by Karina Zikan to graph individual ICESat-2
%%% transects
%%%
%%% SPECIFIED INPUTS:
%%%     icesat2 = path to ICESat-2 file (output from batch_icesat2_terrain_parameter_comparison.m)
%%%     ref_elevations = path to reference elevations file (output from DTM_reference_elevations_calculation.m)
%%%     abbrev = site abriviation for file name
%%% OUTPUTS:
%%%
%%%
%%% Last updated: May 2023 by Karina Zikan

%% Inputs
clearvars;
addpath(['/Users/karinazikan/Documents/ICESat2-AlpineSnow/functions'])
addpath(['/Users/karinazikan/Documents/cmocean'])

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/RCEW/';
%site abbreviation for file names
abbrev = 'RCEW';

% DEM path
DTM_name = [folderpath 'DEMs/RCEW_1m_WGS84UTM11_WGS84.tif'];

%%
%File paths
icesat2_atl08 = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL08-params'];
ref_elevations_atl08 = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL08-ref-elevations'];

icesat2_atl06 = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr'];
ref_elevations_atl06 = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-params'];

icesat2_atl06_class = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06-atl08class'];
ref_elevations_atl06_class = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-atl08class-params'];

%set colors
colors{1} = cmocean('-dense',6);
colors{2} = cmocean('-algae',5);
colors{3} = cmocean('ice',5);

%% Load data ({1} = ATL08, {2} = ATL06, {3} = ATL06 w/ ATL08 classification)
%load the reference elevation data
E{1} = readtable(ref_elevations_atl08);
E{2} = readtable(ref_elevations_atl06);
E{3} = readtable(ref_elevations_atl06_class);

%load the ICESat-2 data
I{1} = readtable(icesat2_atl08); %read in files
I{2} = readtable(icesat2_atl06);
I{3} = readtable(icesat2_atl06_class);

footwidth = 11; % approx. width of icesat2 shot footprint in meters

%DEM
[DTM,Ref] = readgeoraster(DTM_name);
if isfield(Ref,'LatitudeLimits')
    [latgrid,longrid] = meshgrid(Ref.LongitudeLimits(1)+0.5*Ref.CellExtentInLongitude:Ref.CellExtentInLongitude:Ref.LongitudeLimits(2)-0.5*Ref.CellExtentInLongitude,...
        Ref.LatitudeLimits(2)-0.5*Ref.CellExtentInLatitude:-Ref.CellExtentInLatitude:Ref.LatitudeLimits(1)+0.5*Ref.CellExtentInLatitude);
    [xgrid, ygrid,~] = wgs2utm(latgrid,longrid);
else
    x = Ref.XWorldLimits(1)+0.5*Ref.CellExtentInWorldX:Ref.CellExtentInWorldX:Ref.XWorldLimits(2)-0.5*Ref.CellExtentInWorldX;
    if strcmp(Ref.ColumnsStartFrom,'north')
        y = Ref.YWorldLimits(2)-0.5*Ref.CellExtentInWorldY:-Ref.CellExtentInWorldY:Ref.YWorldLimits(1)+0.5*Ref.CellExtentInWorldY;
    else
        y = Ref.YWorldLimits(1)+0.5*Ref.CellExtentInWorldY:Ref.CellExtentInWorldY:Ref.YWorldLimits(2)-0.5*Ref.CellExtentInWorldY;
    end
    [xgrid, ygrid] = meshgrid(x, y); % create grids of each of the x and y coords
end
DTM(DTM>(2.5*10^3)) = NaN;
x = x./(10^3); % to km
y = y./(10^3);% to km
%y = flip(y);
DTM(DTM<0) = NaN;

%% Makesnow-off arrays ({1} = ATL08, {2} = ATL06, {3} = ATL06 w/ ATL08 classification)

%ATL08
% season (1=winter(Jan,Feb,Mar),2=spring(Apr,May,Jun),3=summer(Jul,Aug,Sep),4=fall(Oct,Nov,Dec))
ix_off = find(I{1}.season == 3);
I_off{1} = I{1}(ix_off,:);
E_off{1} = E{1}(ix_off,:);

%ATl06
ix_off = find(I{2}.time.Month >=6 & I{2}.time.Month <= 9);
I_off{2} = I{2}(ix_off,:);
E_off{2} = E{2}(ix_off,:);

%ATl06 atl08class
ix_off = find(I{3}.time.Month >=6 & I{3}.time.Month <= 9);
I_off{3} = I{3}(ix_off,:);
E_off{3} = E{3}(ix_off,:);

%% Calculate Snow-on and Snow-off residuals & snow depth (snow-on residuals shifted by meadian Snow-off residual)
for i = 1:3
    if i == 1
        %         Residuals_off{i} =  I_off{i}.Elevation_bestfit - E_off{i}.elevation_report_fitted;
        Residuals_off{i} =  I_off{i}.Elevation - E_off{i}.elevation_report_fitted;
    else
        Residuals_off{i} =  I_off{i}.h_mean - E_off{i}.elevation_report_nw_mean;
    end
    Residuals_off{i}(Residuals_off{i} > 80) = NaN; Residuals_off{i}(Residuals_off{i} < -80) = NaN; %remove extreme outliers
end
%% Group data
clear SnowDepthTable

n = 3; % ATL06 classified
Residuals =  I{n}.h_mean - E{n}.elevation_report_nw_mean;
Residuals(Residuals > 80) = NaN; Residuals(Residuals < -80) = NaN; %remove extreme outliers
SnowDepthAll = Residuals - median(Residuals_off{n},'omitnan');

E{n}.elevation_report_nw_mean(isnan(SnowDepthAll)) = NaN;

% group elevation
figure(1);
h = histogram(E{n}.elevation_report_nw_mean,30);
elev_binwidth = h.BinWidth; elev_binedges = h.BinEdges;
xlabel('Elevation (m)')
clear h;
for i = 2:length(elev_binedges)
    bins{i-1} = num2str(elev_binedges(i));
end
groupElev = discretize(E{n}.elevation_report_nw_mean, elev_binedges,'categorical',bins);

% group aspect
figure(2);
h = histogram(E{n}.aspect_mean,30);
aspect_binwidth = h.BinWidth; aspect_binedges = h.BinEdges;
xlabel('Aspect (degree)')
clear h;
for i = 2:length(aspect_binedges)
    bins{i-1} = num2str(aspect_binedges(i));
end
groupAspect = discretize(E{n}.aspect_mean, aspect_binedges,'categorical',bins);

% group slope
figure(3);
h = histogram(E{n}.slope_mean,30);
slope_binwidth = h.BinWidth; slope_binedges = h.BinEdges;
xlabel('Slope (degree)')
clear h;
for i = 2:length(slope_binedges)
    bins{i-1} = num2str(slope_binedges(i));
end
groupSlope = discretize(E{n}.slope_mean, slope_binedges,'categorical',bins);

% Make table that can be grouped by various parameters
SnowDepthTable(:,1) = table([datetime(I{n}.time.Year,I{n}.time.Month,I{n}.time.Day)]); %Date
SnowDepthTable(:,2) = table(groupElev);  %binned elevations
SnowDepthTable(:,3) = table(groupAspect);  %binned aspects
SnowDepthTable(:,4) = table(groupSlope);  %binned slopes
SnowDepthTable(:,5) = table(SnowDepthAll); %ATL06_class Snow Depth
SnowDepthTable = renamevars(SnowDepthTable,["Var1","Var2","Var3","Var4","Var5"],["Date","Elevation","Aspect","Slope","ATL06_classSnowDepth"]);

% Average snow depth by paramaters and date
SnowDepthMedTable_elev = varfun(@(x)median(x,'omitnan'),SnowDepthTable,'GroupingVariables',{'Date','Elevation'});
SnowDepthMedTable_aspect = varfun(@(x)median(x,'omitnan'),SnowDepthTable,'GroupingVariables',{'Date','Aspect'});
SnowDepthMedTable_slope = varfun(@(x)median(x,'omitnan'),SnowDepthTable,'GroupingVariables',{'Date','Slope'});
SnowDepthMedTable_date = varfun(@(x)median(x,'omitnan'),SnowDepthTable,'InputVariables','ATL06_classSnowDepth','GroupingVariables','Date');
dates = unique(SnowDepthMedTable_elev.Date);

% Make grouped snow depths arrays
for k = 1:length(dates)
    for i = 2:length(elev_binedges)
        % elevation
        ix = SnowDepthMedTable_elev.Date == dates(k); %index by date
        T = SnowDepthMedTable_elev(ix,:);
        if sum(ismember(num2str(elev_binedges(i)),T.Elevation)) ~= 0
            ix = find(T.Elevation == num2str(elev_binedges(i))); %index by elevation
            SnowDepthMedian = T.Fun_ATL06_classSnowDepth(ix,:);
        else
            SnowDepthMedian = NaN;
        end
        SnowDepthMedArray_elev(k,[i-1]) = SnowDepthMedian;

        % aspect
        ix = SnowDepthMedTable_aspect.Date == dates(k); %index by date
        T = SnowDepthMedTable_aspect(ix,:);
        if sum(ismember(num2str(aspect_binedges(i)),T.Aspect)) ~= 0
            ix = find(T.Aspect == num2str(aspect_binedges(i))); %index by elevation
            SnowDepthMedian = T.Fun_ATL06_classSnowDepth(ix,:);
        else
            SnowDepthMedian = NaN;
        end
        SnowDepthMedArray_aspect(k,[i-1]) = SnowDepthMedian;

        % slope
        ix = SnowDepthMedTable_slope.Date == dates(k); %index by date
        T = SnowDepthMedTable_slope(ix,:);
        if sum(ismember(num2str(slope_binedges(i)),T.Slope)) ~= 0
            ix = find(T.Slope == num2str(slope_binedges(i))); %index by elevation
            SnowDepthMedian = T.Fun_ATL06_classSnowDepth(ix,:);
        else
            SnowDepthMedian = NaN;
        end
        SnowDepthMedArray_slope(k,[i-1]) = SnowDepthMedian;
    end
end
%% Track Plots
i = 5;
date = dates(i);

ix = find(SnowDepthTable.Date == date);
I_track = I{n}(ix,:);
E_track = E{n}(ix,:);
gts = unique(I_track.gt);

for k = 1:length(gts)
    ix = find(I_track.gt == gts(k));
    I_beam = I_track(ix,:);
    E_beam = E_track(ix,:);
    SnowDepth_beam = I_beam.h_mean-E_beam.elevation_report_nw_mean;

    figure(3+k); clf
    subplot(2,2,1); hold on
    scatter([I_beam.Easting./(10^3)],E_beam.elevation_report_nw_mean,'filled')
    scatter([I_beam.Easting./(10^3)],[I_beam.h_mean-median(Residuals_off{3}, 'omitnan')],'filled')
    xlabel('Easting (km)'); ylabel('Elevation (m)');
    legend('Snow Free DEM','ICESat-2')
    set(gca,'fontsize',16);
    set(gcf,'position',[50 50 800 400]);

    subplot(2,2,3); hold on
    stem([I_beam.Easting./(10^3)], SnowDepth_beam,'Linewidth', 2)
    xlabel('Easting (km)'); ylabel('Snow Depth (m)');
    legend('Snow Depth')
    set(gca,'fontsize',16);
    set(gcf,'position',[50 50 800 400]);

    subplot(1,2,2); hold on
    imagesc(x,y,DTM)
    daspect([1 1 1])
    colormap(cmocean('grey'))
    scatter([I_beam.Easting./(10^3)],[I_beam.Northing./(10^3)],[],'g','.')
    xlabel('Easting (km)')
    ylabel('Northing (km)')
    set(gca,'fontsize',16);
    set(gca,'Ydir','normal')
end













