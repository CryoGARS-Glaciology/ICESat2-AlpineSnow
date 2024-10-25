%% SPECIFIED INPUTS:
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
%%% Last updated: Sept 2024 by Karina Zikan & Ellyn Enderlin

%% Inputs
clearvars; close all;
addpath('/bsuhome/karinazikan/scratch/') % path to location of reference_elevations & ICESat2_FootprintCorners functions
%set up the grid size
A1 = -8:8;
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

%ICESat-2 product acronym
acronym = 'A6-40'; %for custom ATL06 with ATL08 classification set to A6-20 for 20m, A6-40 for 40m, 'A6-30' for 30m

%Set output name - MAKE SURE FILENAME SUFIX IS CORRECT!!!!!!!!!!!!!!!!!!!
% file name formats: '-ref-elevations-grid-grad-decent'
%                    '-atl08class-ref-elevations-grid-grad-decent'
%                    '-20m-ref-elevations-grid-grad-decent'
%                    '-atl08class-20m-ref-elevations-grid-grad-decent'
filename_sufix = '-ref-elevations-grid-search-Agg';

%% Set output name
outputname = [abbrev,'-ICESat2-',acronym, filename_sufix, '.csv'];
outputname_acc = [abbrev,'-ICESat2-',acronym, filename_sufix, '-acc.csv'];
outputname_dec = [abbrev,'-ICESat2-',acronym, filename_sufix, '-dec.csv'];

%% Read in files
% Set default_length
if acronym == 'ATL08'
    default_length = 100;
elseif acronym == 'ATL06'
    default_length = 40;
elseif acronym == 'A6-40'
    default_length = 40;
elseif acronym == 'A6-20'
    default_length = 20;
elseif acronym == 'A6-30'
    default_length = 30;
else
    error('acronym must be ATL06 or or A6-40 or A6-20 or ATL08')
end
footwidth = 11; % approx. width of icesat2 shot footprint in meters

%days of year
modays_norm = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays_norm = cumsum(modays_norm); cumdays_norm = [0 cumdays_norm(1:11)];
modays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
cumdays_leap = cumsum(modays_leap); cumdays_leap = [0 cumdays_leap(1:11)];

%read in the snow-off reference elevation map
cd(DTM_path);
if contains(DTM_name,'.tif')
    [DTM,Ref] = readgeoraster(DTM_name);
    DEMdate = DTM_date;
elseif contains(DTM_name,'.mat')
    load(DTM_name);
    for i = 1:length(Z)
        DEMdate(i) = Z(i).deciyear;
    end
end
slope = readgeoraster(DTM_slope);
aspect = readgeoraster(DTM_aspect);

%identify the ICESat-2 data csv files
cd(csv_path);

%filter R2erence DTM elevations
elevations = DTM;
elevations(elevations < -10) = nan; % throw out trash data
elevations(elevations > 10000) = nan; % more trash takeout

%load the ICESat-2 data
T = table; %create a table
icesat2 = [csv_path,csv_name]; %compile the file name
file = readtable(icesat2); %read in files
T = [T; file];

%extract data from columns
zmod = T.h_mean(:); % save the median 'model' elevations (icesat-2 elevations)
zstd = T.h_sigma; %save the standard deviation of the icesat-2 elevation estimates
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings
tracks = T.spot(:); % pull out the beam (1-6)

%identify the ends of each transect and flag them so that neighboring
%transects aren't used when constructing footprints (use beam variable & date)
dates = datetime(T.time.Year,T.time.Month,T.time.Day);
[unique_dates,unique_refs] = unique(dates);
end_flag = zeros(size(norths,1),1);
end_flag(unique_refs) = 1; end_flag(unique_refs(unique_refs~=1)-1) = 1; end_flag(end) = 1;

%filter the data to remove any single-footprint tracks & resave the csv
for k = 1:length(unique_dates)
    ix = find(dates == unique_dates(k));
    if length(ix) == 1
        zmod(ix) = []; zstd(ix) = []; easts(ix) = []; norths(ix) = [];
        tracks(ix) = []; dates(ix) = []; end_flag(ix) = [];
        T(ix,:) = [];
    end
end
%resave icesat2
writetable(T,icesat2);
%find edited unique dates
clear unique_*;
[unique_dates,~] = unique(dates);

%% Snow free data
ix_off = find(T.snowcover == 0);
zmod_off = zmod(ix_off,:);
norths_off = norths(ix_off,:);
easts_off = easts(ix_off,:);
tracks_off = tracks(ix_off,:);

%identify dates and if each footprint is at the end of a track
dates_off = dates(ix_off,:);
[unique_dates_off,unique_refs] = unique(dates_off);
end_flag_off = zeros(size(norths(ix_off,:),1),1);
end_flag_off(unique_refs) = 1; end_flag_off(unique_refs(unique_refs~=1)-1) = 1; end_flag_off(end) = 1;

%% Aggrigated accending and decending corregistation

%determine ascending and descending passes 
zmod_acc = []; zmod_dec = []; zmod_acc_off = []; zmod_dec_off = [];
norths_acc_off = []; easts_acc_off = []; end_flag_acc_off = [];
norths_acc = []; easts_acc = []; end_flag_acc = [];
norths_dec_off = []; easts_dec_off = []; end_flag_dec_off = [];
norths_dec = []; easts_dec = []; end_flag_dec = [];
for k = 1:length(unique_dates)
    for j = 1:6
        track_flag = find(tracks(ix)==j);
        track_norths = norths(ix(track_flag)); track_easts = easts(ix(track_flag));
        track_dirs = (track_norths(2:end)-track_norths(1:end-1))./(track_easts(2:end)-track_easts(1:end-1));
        track_dir(j) = nanmedian(track_dirs);
        clear track_flag track_norths track_easts track_dirs;
    end
    sat_dir(k) = nanmedian(track_dir); clear track_dir;
    %identify data for the date
    ix = find(dates == unique_dates(k));
    ix_off = find(dates_off == unique_dates(k));
    % make accending and decending datasets
    if sat_dir == 0
        zmod_acc = [zmod_acc; zmod(ix)];
        norths_acc = [norths_acc; norths(ix)];
        easts_acc = [easts_acc; easts(ix)];
        end_flag_acc = [end_flag_acc; end_flag(ix)];
        zmod_acc_off = [zmod_acc_off; zmod_off(ix_off)];
        norths_acc_off = [easts_acc_off; easts_off(ix_off)];
        easts_acc_off = [end_flag_acc_off; end_flag_off(ix_off)];
        end_flag_acc_off = [end_flag_acc_off; end_flag_off(ix_off)];
    elseif sat_dir == 1
        zmod_dec = [zmod_dec; zmod(ix)];
        norths_dec = [norths_dec; norths(ix)];
        easts_dec = [easts_dec; easts(ix)];
        end_flag_dec = [end_flag_dec; end_flag(ix)];
        zmod_dec_off = [zmod_dec_off; zmod_off(ix_off)];
        norths_dec_off = [easts_dec_off; easts_off(ix_off)];
        easts_dec_off = [end_flag_dec_off; end_flag_off(ix_off)];
        end_flag_dec_off = [end_flag_dec_off; end_flag_off(ix_off)];
    else
    end
end


% accending coregistration
disp('Accending Coregistration:')
% create the grid search function
GridSearchFunc = @(A)reference_elevations(zmod_acc_off, norths_acc_off, easts_acc_off, end_flag_acc_off, default_length, elevations, slope, aspect, Ref, A); %create the handle to call the coregistration function

%run the grid search
tic
for i = 1:length(A1)
    for j = 1:length(A1)
        fprintf('running cell [ %5.2f , %5.2f ] \n', i, j);
        rmad_grid(i,j) = GridSearchFunc([A1(i),A1(j)]);
    end
    writematrix(rmad_grid,[abbrev,'_',acronym,'_rmadGrid_acc.csv'])
end
toc
% print old and new RMAD
tic
fprintf('x-offset = %5.2f m & y-offset = %5.2f m w/ Accending RNMAD = %5.2f m \n',A1(col),A1(row),min(rmad_grid(:)));
fprintf('Old Accending RNMAD = %5.2f \n', rmad_grid(6,6));
toc

% Calculate corregistered reference elevations
[~,E] = reference_elevations(zmod_acc, norths_acc, easts_acc, end_flag_acc, default_length, elevations, slope, aspect, Ref, [A1(row),A1(col)]); %calculate ref elevations with the shift

% save refelevation csv
writetable(E,outputname_acc);

% Decending coregistration
disp('Decending Coregistration:')
% create the grid search function
GridSearchFunc = @(A)reference_elevations(zmod_dec_off, norths_dec_off, easts_dec_off, end_flag_dec_off, default_length, elevations, slope, aspect, Ref, A); %create the handle to call the coregistration function

%run the grid search
tic
for i = 1:length(A1)
    for j = 1:length(A1)
        fprintf('running cell [ %5.2f , %5.2f ] \n', i, j);
        rmad_grid(i,j) = GridSearchFunc([A1(i),A1(j)]);
    end
    writematrix(rmad_grid,[abbrev,'_',acronym,'_rmadGrid_acc.csv'])
end
toc
% print old and new RMAD
tic
fprintf('x-offset = %5.2f m & y-offset = %5.2f m w/ Decending RNMAD = %5.2f m \n',A1(col),A1(row),min(rmad_grid(:)));
fprintf('Old Decending RNMAD = %5.2f \n', rmad_grid(6,6));
toc

% Calculate corregistered reference elevations
[~,E] = reference_elevations(zmod_dec, norths_dec, easts_dec, end_flag_dec, default_length, elevations, slope, aspect, Ref, [A1(row),A1(col)]); %calculate ref elevations with the shift

% save refelevation csv
writetable(E,outputname_dec);
