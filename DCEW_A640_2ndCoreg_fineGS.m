%%% Code to calculate reference elevations for each ICESat-2 segment with
%%% a 8 by 8 m grid search for aggragated ICESat-2 data. This
%%% code should be performed after DTM_ref_elevations_aggCoreg.m is run
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
%%%     Arow = corse Easting coregistration transform calculated from DTM_ref_elevations_aggCoreg.m
%%%     Acol = corse Northing coregistration transform calculated from DTM_ref_elevations_aggCoreg.m
%%% OUTPUTS:
%%%     Reference_Elevations = csv datatable reporting the non-weighted
%%%         mean, std, weighted mean, and fitted refference elevations,
%%%         mean slope, std slope, mean aspect, std aspect
%%%
%%%
%%% Last updated: June 2025 by Karina Zikan

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

% Corse Grid search offset (from corse coregistration)
Arow = 0;
Acol = 0;

%% Set output name
filename_sufix = 'A6-40-ref-elevations-grid-search-fineGS-agg';
outputname = [abbrev,'-ICESat2-', filename_sufix, '.csv'];

%ICESat-2 product acronym
acronym = 'A6-40';
default_length = 40; %segment length
footwidth = 11; %average segment width

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

zmod = T.h_mean(:); % mean elevations (icesat-2 elevations)
zstd = T.h_sigma; %standard deviation of the icesat-2 elevation estimates
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings

%% Snow free data
ix_off = find(T.snowcover == 0);

%identify the ends of each transect and flag them so that neighboring
%transects aren't used when constructing footprints (use beam variable & date)
[~,unique_refs] = unique([convertTo(dates, 'yyyymmdd'), T.gt],'rows');
[unique_dates,~] = unique(dates);
end_flag = zeros(size(norths,1),1);
end_flag(unique_refs) = 1; end_flag(unique_refs(unique_refs~=1)-1) = 1; end_flag(end) = 1;

dates_off = dates(ix_off,:);
gt_off = T.gt(ix_off,:);
[~,unique_refs] = unique([convertTo(dates_off, 'yyyymmdd'), gt_off],'rows');
[unique_dates_off,~] = unique(dates_off);
end_flag_off = zeros(size(norths(ix_off,:),1),1);
end_flag_off(unique_refs) = 1; end_flag_off(unique_refs(unique_refs~=1)-1) = 1; end_flag_off(end) = 1;

%reference elevation function to optimise
GradDecentFunc = @(A)reference_elevations(zmod(ix_off,:), norths(ix_off,:), easts(ix_off,:), end_flag_off, default_length, elevations, slope, aspect, Ref, A); %create the handle to call the coregistration function

%% Grid of possible inputs to calculate initial guess
A1 = -0.9:0.1:0.9;

tic
for i = 1:length(A1)
    tic
    for j = 1:length(A1)
        rmad_grid(i,j) = GradDecentFunc([Arow+A1(i),Acol+A1(j)]);
    end
    fprintf('row %5.2f completed',i);
    writematrix(rmad_grid,[abbrev,'_rmadFineGrid.csv'])
    toc
end
toc


figure(2);
im = imagesc(rmad_grid);
xticks(1:length(A1)); yticks(1:length(A1));
xticklabels(A1); yticklabels(A1);
colorbar;

[row, col] = find(ismember(rmad_grid, min(rmad_grid(:))));
A1row = A1(row); A1col = A1(col);

writematrix([Arow+A1row,Acol+A1col],[abbrev,'_Ashift.csv'])

% Test Gradient Decent from grid guess
tic
fprintf('x-offset = %5.2f m & y-offset = %5.2f m w/ RNMAD = %5.2f m \n',Arow+A1row,Acol+A1col,min(rmad_grid(:)));
fprintf('Old RNMAD = %5.2f', rmad_grid(10,10));
toc

%% Calculate corregistered reference elevations
[~,E] = reference_elevations(zmod, norths, easts, end_flag, default_length, elevations, slope, aspect, Ref, [Arow+A1row,Acol+A1col]); %calculate ref elevations with the shift

%% save ref elevation csv
writetable(E,outputname);