%%% Code to calculate reference elevations for each ICESat-2 segment with
%%% no coregistration transform
%%%
%%% Writen by Karina Zikan & Ellyn Enderlin
%%%
%%% SPECIFIED INPUTS:
%%%     DTM_path = path to the reference DTM on your computer
%%%     DTM_name = DTM file name
%%%     DTM_slope = Slope map file name
%%%     DTM_aspect = Aspect map file name
%%%     csv_path = path to the ICESat-2 datafiles on your computer
%%%     csv_name = name of ICESat-2 csv file
%%%     abbrev = site abriviation for file name
%%% OUTPUTS:
%%%     Reference_Elevations = csv datatable reporting the non-weighted
%%%         mean, std, weighted mean, and fitted refference elevations,
%%%         mean slope, std slope, mean aspect, std aspect
%%%
%%%
%%% Last updated: Jun 2025 by Karina Zikan

%% Inputs
clearvars; close all;
addpath('/bsuhome/karinazikan/scratch/') % path to location of reference_elevations & ICESat2_FootprintCorners functions

%DTM (be sure the path ends in a /)
DTM_path = '/bsuhome/karinazikan/scratch/DCEW/'; %path to dtm, slope, & aspect maps
DTM_name = 'DryCreekBase1m_WGS84UTM11_DEM.tif';
if contains(DTM_name,'.tif')
    DTM_date = '20120826'; %only need to change this if the DTM is a geotiff
end
% Slope
DTM_slope = 'DryCreekBase1m_WGS84UTM11-slope-002.tif';
% Aspect
DTM_aspect = 'DryCreekBase1m_WGS84UTM11-aspect-001.tif';

%csv (be sure the path ends in a /)
csv_path = '/bsuhome/karinazikan/scratch/DCEW/A6-40/'; %Path to ICESat-2 data with snow cover classification
csv_name = 'DCEW-ICESat2-A6-40-SnowCover.csv';

%site abbreviation for file names
abbrev = 'DCEW';


%% Set output name
filename_sufix = 'A6-40-ref-elevations-noCoreg';
outputname = [abbrev,'-ICESat2-', filename_sufix, '.csv'];

%ICESat-2 product acronym
acronym = 'A6-40';
default_length = 40; %segment length
footwidth = 11; %average segment width

% no trasform for no coregistration
Arow = 0;
Acol = 0;

%% Read in files
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

%filter R2erence DTM elevations
elevations = DTM;
elevations(elevations < -10) = nan; % throw out trash data
elevations(elevations > 10000) = nan; % more trash takeout

%load the ICESat-2 data
T = table; %create a table
icesat2 = [csv_path,csv_name]; %compile the file name
file = readtable(icesat2); %read in files
T = [T; file];
%T = T(1:5000,:); % ONLY FOR TESTING!!!!!!!!!!

% sort time and gt in T
T = sortrows(T,{'time','gt'});
dates = datetime(T.time.Year,T.time.Month,T.time.Day);

% pull out needed variables
zmod = T.h_mean(:); %mean elevations (icesat-2 elevations)
zstd = T.h_sigma; %standard deviation of the icesat-2 elevation estimates
easts = T.Easting(:); % easting values
norths = T.Northing(:); % northing values

%% Snow free data
ix_off = find(T.snowcover == 0);

% Identify the ends of each transect (using beam variable & date) and flag them so that neighboring
% transects aren't used when constructing footprints
[~,unique_refs] = unique([convertTo(dates, 'yyyymmdd'), T.gt],'rows');
end_flag = zeros(size(norths,1),1);
end_flag(unique_refs) = 1; end_flag(unique_refs(unique_refs~=1)-1) = 1; end_flag(end) = 1;

%% Calculate corregistered reference elevations
[~,E] = reference_elevations(zmod, norths, easts, end_flag, default_length, elevations, slope, aspect, Ref, [Arow,Acol]); %calculate ref elevations

%% save ref elevation csv
writetable(E,outputname);
