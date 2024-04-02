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
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/MCS/';
%site abbreviation for file names
abbrev = 'MCS';

% DEM path
DTM_name = [folderpath 'DEMs/MCS_REFDEM_WGS84.tif'];
% DTM_name = [folderpath 'DEMs/RCEW_1m_WGS84UTM11_WGS84.tif'];
% DTM_name = [folderpath 'DEMs/Banner_Bare_Earth_DEMs_mosaic_UTM11WGS84.tif'];
% DTM_name = [folderpath 'DEMs/DryCreekBase1m_WGS84UTM11_DEM.tif'];

%%
%File paths
icesat2_atl06_class = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-atl08class'];
ref_elevations_atl06_class = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-atl08class-params'];

%% Load data ({1} = ATL08, {2} = ATL06, {3} = ATL06 w/ ATL08 classification)
%load the reference elevation data
filepath = strcat(folderpath, 'IS2_Data/ATL06-atl08class-AllData.csv');
df = readtable(filepath);
df.Easting = df.Easting./(10^3);
df.Northing = df.Northing./(10^3);

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


%% Track Plots
date = '23-Feb-2020';
dates = datetime(df.time.Year,df.time.Month,df.time.Day);
ix = find(dates == date);
Track = df(ix,:);

figure(1); clf
imagesc(x,y,DTM); hold on
daspect([1 1 1]); colormap(cmocean('grey'));
scatter([df.Easting],[df.Northing],[],'green','.')
scatter([Track.Easting./(10^3)],[Track.Northing./(10^3)],[],'magenta','.')
xlabel('Easting [km]'); ylabel('Northing [km]');
set(gca,'fontsize',20); set(gca,'Ydir','normal');
c = colorbar;
c.Label.String = 'Elevation (m)';

% %% Residual Map Plots
% figure(2); clf
% imagesc(x,y,DTM); hold on
% daspect([1 1 1]); colormap(cmocean('grey'));
% scatter(df,'Easting','Northing','elev_residuals_vertcoreg_slopecorrected')
% xlabel('Easting [km]'); ylabel('Northing [km]');
% set(gca,'fontsize',20); set(gca,'Ydir','normal');
% c = colorbar;
% c.Label.String = 'Elevation (m)';
% 












