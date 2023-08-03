%%% This code callculates a 11m raster from ICESat-2 ATL03 points by taking
%%% the minimum ATL03 point within each 11m cell. Geolocation information
%%% comes from the the site reference DTM
%%%
%%% SPECIFIED INPUTS:
%%%     folderpath = path to the site folder on your computer
%%%     DTM_name = DTM file name
%%%     abbrev = site abriviation for file name
%%% OUTPUTS:
%%%     IS2grid = 11m geotif raster
%%% 
%%% Last updated: Aug 2023 by Karina Zikan
clear all

%% Inputs
%Folder path 
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/RCEW/';
DTM_name = 'RCEW_1m_WGS84UTM11_WGS84.tif';
%site abbreviation for file names
abbrev = 'RCEW';

%% Read in files
%Read in files
file = readtable([folderpath 'IS2_Data/' abbrev '-ICESat2-ATL03']);
[DTM,Ref] = readgeoraster([folderpath 'DEMs/' DTM_name]);

% crs check
if Ref.ProjectedCRS.Name ~= 'WGS 84 / UTM zone 11N'
    error('Error: DTM must be in WGS84 / UTM zone 11N')
end

%DEMgrid
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


%% visualize data
figure(1);
imagesc(x,y,DTM);
hold on
scatter3(file,"X","Y","Z",Marker=".");

%% Make 11m elevation array
x11 = downsample(x,11); y11 = downsample(y,11);

for i = 1:length(x11)
    for j = 1:length(y11)
        ix = find(file.X <= (x11(i)+ 5.5 & file.X >= (x11(i)-5.5) & ...
            file.Y <= (y11(i)+ 5.5 & file.Y >= (y11(i)-5.5)))); % Identify ATL03 points w/in cell
        if isempty(ix) == 1
            IS2grid(j,i) = NaN; %Set elevation to NaN for empty cells
        else
            IS2grid(j,i) = min(file.Z(ix)); %Set elevation to min ATL03
        end
    end
end

%% Create the referencing matrix (R)
R = maprasterref('RasterSize',size(IS2grid), ...
    'XLimWorld',[min(x11) max(x11)],'YLimWorld',[min(y11) max(y11)]);
R.ProjectedCRS  = Ref.ProjectedCRS; 
R.ColumnsStartFrom = 'north';

%% Export raster
filename = [folderpath 'DEMs/' abbrev '-ICESat2-ATL03.tif'];
geotiffwrite(filename,IS2grid,R,'CoordRefSysCode','epsg:32611')