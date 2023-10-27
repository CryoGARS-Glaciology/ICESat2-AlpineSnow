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
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/Banner/';
DTM_name = 'Banner_Bare_Earth_DEMs_mosaic_UTM11WGS84.tif';
%site abbreviation for file names
abbrev = 'Banner';

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
scatter3(file.X,file.Y,file.Z,Marker=".");


%% Make 11m elevation array

x11 = downsample(x,11); y11 = downsample(y,11);

for i = 1:length(x11)
    for j = 1:length(y11)
        ix = find(file.X <= (x11(i)+ 5.5) & file.X >= (x11(i)-5.5) & file.Y <= (y11(j)+ 5.5) & file.Y >= (y11(j)- 5.5));
        check(j,i) = sum(ix);
        if isempty(ix) == 1
            IS2grid(j,i) = NaN; %Set elevation to NaN for empty cells
        else

            #IS2grid(j,i) = nanmean(file.Z(ix));
            IS2grid(j,i) = min(file.Z(ix)); %Set elevation to min ATL03

        end
    end
end

%% Make 11m ref elevation array
for i = 1:length(x11)
    for j = 1:length(y11)
        ix = find(x <= (x11(i)+5.5) & x >= (x11(i)-5.5));
        iy = find(y <= (y11(j)+5.5) & y >= (y11(j)-5.5));
        DTM_11(j,i) = nanmean(DTM(iy,ix),'all');
    end
end
DTM_11 = double(DTM_11);

%% Create the referencing matrix (R)
R = maprasterref('RasterSize',size(IS2grid), ...
    'XLimWorld',[min(x11) max(x11)],'YLimWorld',[min(y11) max(y11)]);
R.ProjectedCRS  = Ref.ProjectedCRS; 
R.ColumnsStartFrom = 'north';

%% Export ATL03 raster
filename = [folderpath 'DEMs/' abbrev '-ICESat2-ATL03-test.tif'];
geotiffwrite(filename,IS2grid,R,'CoordRefSysCode','epsg:32611')

%% Export 11m DTM raster
filename = [folderpath 'DEMs/RCEW_1m_WGS84UTM11_WGS84_11m.tif'];
geotiffwrite(filename,DTM_11,R,'CoordRefSysCode','epsg:32611')

%% Save IS2grid as csv just in case
writematrix(IS2grid,"RCEW_ATL03_grid_mean.csv");

%% Save DTM_11 as csv just in case
writematrix(DTM_11,"DTM_11.csv");
%%
for i = 1:length(ix)
    for j = 1:length(ix)
        check(i,j) = sum(ix{i,j});
    end 
end
