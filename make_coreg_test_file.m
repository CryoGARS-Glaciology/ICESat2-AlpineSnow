clear all

%% inputs
csv_path = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/RCEW/IS2_Data/';
csv_name = 'RCEW-ICESat2-ATL06sr-ref-elevations';

%% load the file
file = readtable([csv_path,csv_name]); %read in files

%% Make shifted test datafile
elevation = file.elevation_report_nw_mean;
