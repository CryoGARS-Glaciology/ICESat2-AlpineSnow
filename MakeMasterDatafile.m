%%% This code was writen by Karina Zikan to calculate the elevation
%%% residuals and compile all the data into one master array
%%%
%%% SPECIFIED INPUTS:
%%%     icesat2 = path to ICESat-2 file (output from batch_icesat2_terrain_parameter_comparison.m)
%%%     ref_elevations = path to reference elevations file (output from DTM_reference_elevations_calculation.m)
%%%     abbrev = site abriviation for file name
%%%     outputname = output file name
%%%
%%% OUTPUTS:
%%%     File = data file with all the data and calcualted residuals
%%%
%%% Last updated: Nov 2023 by Karina Zikan

%% Inputs
clearvars;

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/MCS/';
%site abbreviation for file names
abbrev = 'MCS';

%File paths
icesat2 = [folderpath 'IS2_Data/A6-40/' abbrev '-ICESat2-A6-40-SnowCover'];
ref_elevations = [folderpath 'IS2_Data/A6-40/' abbrev '-ICESat2-A6-40-ref-elevations-noCoreg'];

%output file name
outputname = 'ATL06-A6-40-AllData-noCoreg';

%set footprint length
% footprint = 40;

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

%Calculate ICESat-2 along track slope
IS2_slope = abs(atand(I.dh_fit_dx));

%calculate quadratic slope correction w/ dtm slope
slope = E_off.slope_mean;
x= slope; y = Residuals(ix);
ind = isnan(x) | isnan(y); %index nans
x(ind) = []; y(ind) = []; %remove nans
p = polyfit(x,y,2); % fit quadratic
% Vertical corregistration
elev_residuals_dtm_slopecorrected = Residuals-polyval(p,E.slope_mean);
elev_residuals_dtm_slopecorrected_off = elev_residuals_dtm_slopecorrected(ix);

%calculate quadratic slope correction w/ IS2 slope
slope = IS2_slope(ix);
x= slope; y = Residuals(ix);
ind = isnan(x) | isnan(y); %index nans
x(ind) = []; y(ind) = []; %remove nans
p = polyfit(x,y,2); % fit quadratic
% Vertical corregistration
elev_residuals_is2_slopecorrected = Residuals-polyval(p,IS2_slope);
elev_residuals_is2_slopecorrected_off = elev_residuals_is2_slopecorrected(ix);

% shift by median snow-off residual to vertically coregister residuals
Residuals_VertCoreg = Residuals - median(Residuals_off,'omitnan');
elev_residuals_vertcoreg_dtm_slopecorrected = elev_residuals_dtm_slopecorrected - median(elev_residuals_dtm_slopecorrected_off,'omitnan');
elev_residuals_vertcoreg_is2_slopecorrected = elev_residuals_is2_slopecorrected - median(elev_residuals_is2_slopecorrected_off,'omitnan');

%% Combine all data into one file
Output = [I E];
Output.elev_residuals = Residuals;
Output.elev_residuals_vertcoreg = Residuals_VertCoreg;
Output.elev_residuals_vertcoreg_dtm_slopecorrected = elev_residuals_vertcoreg_dtm_slopecorrected;
Output.elev_residuals_vertcoreg_is2_slopecorrected = elev_residuals_vertcoreg_is2_slopecorrected;
Output.IS2_slope_deg = IS2_slope;

Output_off = Output(Output.snowcover == 0,:);
Output_on = Output(Output.snowcover == 1,:);

%% Write output file
writetable(Output,[folderpath 'IS2_Data/A6-40/' outputname '.csv']);
writetable(Output_off,[folderpath 'IS2_Data/A6-40/' outputname '_off.csv']);
writetable(Output_on,[folderpath 'IS2_Data/A6-40/' outputname '_on.csv']);
