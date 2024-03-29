%%% This code was writen by Karina Zikan and Ellyn Enderlin to calculate
%%% non-weighted mean, weighted mean, and weighted fitted refference
%%% elevations from a reference DTM and mean slope and aspect from slope 
%%% and aspect maps to compare to ICESat-2 ATL08 data
%%% (Use DTM_reference_elevations_calculation_sliderulecsv.m for ATL06 data from sliderule)
%%% << Adding mean vegitation is planned >>
%%% 
%%% SPECIFIED INPUTS:
%%%     DTM_path = path to the reference DTM on your computer
%%%     DTM_name = DTM file name
%%%     DTM_slope = Slope map file name
%%%     DTM_aspect = Aspect map file name
%%%     csv_path = path to the ICESat-2 datafiles on your computer
%%%     csv_name = name of ICESat-2 csv file
%%%     abbrev = site abriviation for file name
%%%     acronym = ICESat-2 product acronym
%%% OUTPUTS:
%%%     Reference_Elevations = csv datatable reporting the non-weighted
%%%         mean, std, weighted mean, and fitted refference elevations,
%%%         mean slope, std slope, mean aspect, std aspect
%%%         
%%%
%%% Last updated: feb 2024 by Karina Zikan


%% Inputs
clearvars; close all;
addpath('./functions') 

%DTM (be sure the path ends in a /)
DTM_path = 'Sites/RCEW/DEMs/';

DTM_name = 'RCEW_1m_WGS84UTM11_WGS84_CoReg.tif';

if contains(DTM_name,'.tif')
    DTM_date = '20120826'; %only need to change this if the DTM is a geotiff
end
% Slope
DTM_slope = 'RCEW_1m_WGS84UTM11_WGS84-slope_CoReg.tif';
% Aspect
DTM_aspect = 'RCEW_1m_WGS84UTM11_WGS84-aspect_CoReg.tif';



%csv (be sure the path ends in a /)
csv_path = '/Users/karinazikan/Documents/GitHub/ICESat2-AlpineSnow/Sites/RCEW/IS2_Data/';
csv_name = 'RCEW-ICESat2-ATL08-params.csv';

%site abbreviation for file names
abbrev = 'RCEW';

%ICESat-2 product acronym
acronym = 'ATL06'; %set to ATL06-20 for the 20m atl06 data

%Set output name - MAKE SURE FILENAME SUFIX IS CORRECT!!!!!!!!!!!!!!!!!!!
filename_sufix = '-ref-elevations-mean-CoReg';

%% Set output name
outputname = [abbrev,'-ICESat2-',acronym, filename_sufix, '.csv'];

%% Read in files
% Set default_length
if acronym == 'ATL08'
    default_length = 100;
elseif acronym == 'ATL06'
    default_length = 40;
elseif acronym == 'ATL06-20'
    default_length = 20;
else
    error('acronym must be ATL06 or ATL06-20 or ATL08')
end 

%days of year
modays_norm = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays_norm = cumsum(modays_norm); cumdays_norm = [0 cumdays_norm(1:11)];
modays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
cumdays_leap = cumsum(modays_leap); cumdays_leap = [0 cumdays_leap(1:11)];

%read in the snow-off reference elevation map
cd_to_DTM = ['cd ',DTM_path]; eval(cd_to_DTM);
if contains(DTM_name,'.tif')
    [DTM,Ref] = readgeoraster(DTM_name);
    DEMdate = DTM_date;
elseif contains(DTM_name,'.mat')
    load_DTM = ['load ',DTM_name]; eval(load_DTM);
    for i = 1:length(Z)
        DEMdate(i) = Z(i).deciyear;
    end
end
slope = readgeoraster(DTM_slope);
aspect = readgeoraster(DTM_aspect);

%identify the ICESat-2 data csv files
cd_to_csv = ['cd ',csv_path]; eval(cd_to_csv);
csvs = dir([acronym,'*.csv']); %if running more than once, rename original csvs to '*raw.csv' and change the search here to find files with that ending

%filter R2erence DTM elevations
% elevations(elevations < -10) = nan;
elevations = DTM;
elevations(elevations < -10) = nan; % throw out trash data
elevations(elevations > 10000) = nan; % more trash takeout

%load the ICESat-2 data
T = table; %create a table
icesat2 = [csv_path abbrev '-ICESat2-ATL08-params']; %compile the file name
file = readtable(icesat2); %read in files
T = [T; file];
%T = T(1:5,:); % ONLY FOR TESTING!!!!!!!!!!

% T = T([1:250],:);
zmod = T.Elevation(:); % save the median 'model' elevations (icesat-2 elevations)
zmodfit = T.Elevation_bestfit(:); % save the fitted 'model' elevations (icesat-2 elevations_bestfit)
zmodfit(isnan(zmod)) = NaN;
zstd = T.std; %save the standard deviation of the icesat-2 elevation estimates
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings
footwidth = 11; % approx. width of icesat2 shot footprint in meters

%identify the ends of each transect and flag them so that neighboring
%transects aren't used when constructing footprints (use beam variable & date)
dates = T.date;
[~,unique_refs] = unique([num2str(dates)],'rows');
end_flag = zeros(size(norths,1),1);
end_flag(unique_refs) = 1; end_flag(unique_refs(unique_refs~=1)-1) = 1; end_flag(end) = 1;


%% Calculating footprints for each data point
%define the Reference elevation data
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

% calculates footprint corners
[xc,yc,theta] = ICESat2_FootprintCorners(norths,easts,default_length,end_flag);

%% Calculate Reference Elevations, Slope, & Aspect
tic
for r=1:length(zmod)

    %identify the R2erence elevation points in each ICESat2 footprint
    xv = xc(r,[3:6 3]); % bounding box x vector
    yv = yc(r,[3:6 3]); % bounding box y vector

    % subset giant grid
    ix = find(x <= (xc(r,1)+100) & x >= (xc(r,1)-100)); % x index for subgrid
    iy = find(y <= (yc(r,1)+100) & y >= (xc(r,1)-100)); % y index for subgrid
    xsubgrid = xgrid(iy,ix);
    ysubgrid = ygrid(iy,ix);
    subelevations = elevations(iy,ix);
    subslope = slope(iy,ix);
    subaspect = aspect(iy,ix);

   %data in the footprint
    in = inpolygon(xsubgrid, ysubgrid, xv, yv); % get logical array of in values
    pointsinx = xsubgrid(in); % save x locations
    pointsiny = ysubgrid(in); % save y locations
    elevationsin = subelevations(in); % save elevations
    slopesin = subslope(in); % save slopes
    aspectsin = subaspect(in); % save slopes

    %wieghted average
    dist = nan([1,length(pointsinx)])'; %initialize dist
    for a = 1:length(pointsinx)
        phi = atan2d((pointsiny(a)-norths(r)),(pointsinx(a)-easts(r)));
        dist(a)=abs(sqrt((pointsiny(a)-norths(r))^2+(pointsinx(a)-easts(r))^2)*sind(phi-theta(r))); %distance from the line in the center of the window
    end
    maxdist = footwidth/2; % defining the maximum distance a point can be from the center icesat2 point
    w = 15/16*(1-(dist/maxdist).^2).^2; %bisqared kernel
    elevation_report_mean(r,:) = sum(w.*elevationsin)./sum(w); %weighted elevation estimate
    elevation_report_std(r,:) = std(elevationsin); %std of the elevations within the footprint

    %non wieghted average
    elevation_report_nw_mean(r,:) = nanmean(elevationsin); % non-wieghted elevations
    slope_mean(r,:) = nanmean(slopesin);
    slope_std(r,:) = std(slopesin);
    aspect_mean(r,:) = nanmean(aspectsin);
    aspect_std(r,:) = std(aspectsin);
end
%interpolated elevation
elevation_report_interp = interp2(x,y,elevations,easts,norths);
toc

% Write reference elevation table
E = table(elevation_report_nw_mean,elevation_report_mean,elevation_report_interp,elevation_report_std,slope_mean,slope_std,aspect_mean,aspect_std);

%% write table to csv
writetable(E,outputname);


%% Sanity Checks
% % Distance and weighting check
% figure;
% plot3(pointsinx, pointsiny,dist,'.')
% xlabel('Easting (km)')
% ylabel('Northing (km)')
% zlabel('Distance from ICESat-2 track centerline (m)')
%
% figure;
% plot3(pointsinx, pointsiny,w,'.')
% xlabel('Easting (km)')
% ylabel('Northing (km)')
% zlabel('Point weighting (m)')




