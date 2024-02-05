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
addpath(['./functions'])
addpath(['/Users/karinazikan/Documents/cmocean'])

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/RCEW/';
%site abbreviation for file names
abbrev = 'RCEW';

% DEM path
DTM_name = [folderpath 'DEMs/RCEW_1m_WGS84UTM11_WGS84.tif'];

%%
%File paths
icesat2_atl06_class = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-atl08class'];
ref_elevations_atl06_class = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-atl08class-params'];

%% Load data ({1} = ATL08, {2} = ATL06, {3} = ATL06 w/ ATL08 classification)
%load the reference elevation data
E = readtable(ref_elevations_atl06_class);

%load the ICESat-2 data
I = readtable(icesat2_atl06_class);

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
dates = datetime(I.time.Year,I.time.Month,I.time.Day);
ix = find(dates == date);
I_track = I(ix,:);
E_track = E(ix,:);


figure(1); clf
imagesc(x,y,DTM); hold on
daspect([1 1 1]); colormap(cmocean('grey'));
scatter([I.Easting./(10^3)],[I.Northing./(10^3)],[],'green','.')
scatter([I_track.Easting./(10^3)],[I_track.Northing./(10^3)],[],'magenta','.')
xlabel('Easting [km]'); ylabel('Northing [km]');
set(gca,'fontsize',20); set(gca,'Ydir','normal');
c = colorbar;
c.Label.String = 'Elevation (m)';













