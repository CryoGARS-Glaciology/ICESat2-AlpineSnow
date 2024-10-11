%% Inputs
clearvars;

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/RCEW/';
%site abbreviation for file names
abbrev = 'RCEW';

%File paths
icesat2 = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06-atl08class-SnowCover'];
ref_elevations = [folderpath 'IS2_Data/' abbrev '-ICESat2-ATL06sr-atl08class-ref-elevations'];


%% Load data
%load the reference elevation data
E = readtable(ref_elevations);

%load the ICESat-2 data
I = readtable(icesat2);

%% Make snow-off array
ix = find(I.snowcover == 0);
I_off = I(ix,:);
E_off = E(ix,:);


%% Calculate residuals w/ and w/out vertical cogregistration
Residuals_off =  I_off.h_mean - E_off.elevation_report_nw_mean;
Residuals =  I.h_mean - E.elevation_report_nw_mean;

%remove extreme outliers
Residuals_off(Residuals_off > 80) = NaN; Residuals_off(Residuals_off < -80) = NaN;
Residuals(Residuals > 80) = NaN; Residuals(Residuals < -80) = NaN;

nmad = 1.4826*median(abs(Residuals_off-nanmean(Residuals_off)),'omitnan')
    %below slope cutoff

