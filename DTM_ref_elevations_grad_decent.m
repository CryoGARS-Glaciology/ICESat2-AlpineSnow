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
%%% Last updated: March 2024 by Karina Zikan


%% Inputs
clearvars; close all;
addpath('./functions') 

%DTM (be sure the path ends in a /)
DTM_path = 'Sites/RCEW/DEMs/';

DTM_name = 'RCEW_1m_WGS84UTM11_WGS84.tif';

if contains(DTM_name,'.tif')
    DTM_date = '20120826'; %only need to change this if the DTM is a geotiff
end
% Slope
DTM_slope = 'RCEW_1m_WGS84UTM11_WGS84-slope.tif';
% Aspect
DTM_aspect = 'RCEW_1m_WGS84UTM11_WGS84-aspect.tif';



%csv (be sure the path ends in a /)
csv_path = '/Users/karinazikan/Documents/GitHub/ICESat2-AlpineSnow/Sites/RCEW/IS2_Data/';
csv_name = 'RCEW-ICESat2-ATL06-atl08class-SnowCover.csv';

%site abbreviation for file names
abbrev = 'RCEW';

%ICESat-2 product acronym
acronym = 'ATL06'; %set to ATL06-20 for the 20m atl06 data

%Set output name - MAKE SURE FILENAME SUFIX IS CORRECT!!!!!!!!!!!!!!!!!!!
    % file name formats: '-ref-elevations-grid-grad-decent' 
    % '-atl08class-ref-elevations-grid-grad-decent' 
    % '-20m-ref-elevations-grid-grad-decent' 
    % '-atl08class-20m-ref-elevations-grid-grad-decent' 
filename_sufix = '-atl08class-ref-elevations-grid-grad-decent';

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


%% ICESat-2 variables
if acronym == 'ATL06'| acronym == 'ATL06-20'
    zmod = T.h_mean(:); % save the median 'model' elevations (icesat-2 elevations)
    zstd = T.h_sigma; %save the standard deviation of the icesat-2 elevation estimates
elseif acronym == 'ATL08'
    zmod = T.Elevation(:); % save the median 'model' elevations (icesat-2 elevations)
    zstd = T.std; %save the standard deviation of the icesat-2 elevation estimates
else
    error('Error: the variable acronym must be ATL06, ATL06-20, or ATL08. Set acronym in the inputs')
end
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings
footwidth = 11; % approx. width of icesat2 shot footprint in meters

%% Snow free data
ix_off = find(T.snowcover == 0);

%identify the ends of each transect and flag them so that neighboring
%transects aren't used when constructing footprints (use beam variable & date)
dates = datetime(T.time.Year,T.time.Month,T.time.Day);
[~,unique_refs] = unique(dates);
end_flag = zeros(size(norths,1),1);
end_flag(unique_refs) = 1; end_flag(unique_refs(unique_refs~=1)-1) = 1; end_flag(end) = 1;

dates_off = dates(ix_off,:);
[~,unique_refs] = unique(dates_off);
end_flag_off = zeros(size(norths(ix_off,:),1),1);
end_flag_off(unique_refs) = 1; end_flag_off(unique_refs(unique_refs~=1)-1) = 1; end_flag_off(end) = 1;
% 
% %% Gradient Decent 
% tic
GradDecentFunc = @(A)reference_elevations(zmod(ix_off,:), norths(ix_off,:), easts(ix_off,:), end_flag_off, default_length, elevations, slope, aspect, Ref, A); %create the handle to call the coregistration function
% [Abest,RMADbest] = fminsearch(GradDecentFunc,[0,0]); %initial horizontal offset estimate = [0,0] = [0 m East, 0 m North]
% %save([abbrev,'-Abest.mat'],'Abest','-v7.3'); save([abbrev,'-RMADbest.mat'],'RMADbest','-v7.3'); %save offset for future
% fprintf('x-offset = %5.2f m & y-offset = %5.2f m w/ RNMAD = %5.2f m \n',Abest(:,1),Abest(:,2),RMADbest); 
% toc

%% Grid of posible inputs to calculate initial guess
A1 = -5:5;

tic 
for i = 1:length(A1)
    tic
    parfor j = 1:length(A1)
        rmad_grid(i,j) = GradDecentFunc([A1(i),A1(j)]);
    end
    toc
    writematrix(rmad_grid,[abbrev,'_rmadGrid.csv'])
end
toc

%%
tic 
for i = 1:length(A1)
    for j = 1:length(A1)
        rmad_grid2(i,j) = GradDecentFunc([A1(i),A1(j)]);
    end
end
toc

figure(2);
imagesc(rmad_grid); 

xticks(1:41); yticks(1:41); 
xticklabels(A1); yticklabels(A1); 
colorbar;

[row, col] = find(ismember(rmad_grid, min(rmad_grid(:))));
Arow = A1(row); Acol = A1(col);
writematrix([Arow,Acol],[abbrev,'_Ashift.csv'])

% %% Test Gradient Decent from grid guess
% tic
% GradDecentFunc = @(A)reference_elevations(zmod(ix_off,:), norths(ix_off,:), easts(ix_off,:), end_flag_off, default_length, elevations, slope, aspect, Ref, A); %create the handle to call the coregistration function
% [Agrid_best,RMADgrid_best] = fminsearch(GradDecentFunc,[Arow,Acol],optimset('PlotFcns',@optimplotfval)); %initial horizontal offset estimate = [0,0] = [0 m East, 0 m North]
fprintf('x-offset = %5.2f m & y-offset = %5.2f m w/ RNMAD = %5.2f m \n',Acol,Arow,min(rmad_grid(:)));
% toc


%% Calculate corregistered reference elevations
[~,E] = reference_elevations(zmod, norths, easts, end_flag, default_length, elevations, slope, aspect, Ref, [Arow,Acol]); %calculate ref elevations with the shift

%% save refelevation csv
writetable(E,outputname);








