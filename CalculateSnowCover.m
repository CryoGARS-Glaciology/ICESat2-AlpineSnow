%%% This code was writen by Karina Zikan to pull snow-on snow-free
%%% information from classified NDSI images
%%%
%%% SPECIFIED INPUTS:
%%%     folderpath = path to site folder
%%%     IS2file = icesat2 csv filename
%%%     abbrev = site abreviation
%%%     acronym = ICESat-2 product acronym (ATL08, ALT06, ect)
%%%     outputname = output file name
%%%
%%% OUTPUTS:
%%%     
%%%
%%% Last updated: Nov 2023 by Karina Zikan

%% Inputs
clearvars;

%Folder path 
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/RCEW/';
SnowCoverPath = [folderpath 'SnowCover/'];
%site abbreviation for file names
abbrev = 'RCEW';

%filename
IS2file = 'RCEW-ICESat2-ATL06sr.csv';
acronym = 'ATL06';

%Set output name
outputname = [abbrev,'-ICESat2-',acronym,'-SnowCover.csv'];

%% read in ICESat-2 data
icesat2path = [folderpath 'IS2_Data/' IS2file];
IS2 = readtable(icesat2path);

%% Set up date arrays
IS2dates = datetime(IS2.time.Year,IS2.time.Month,IS2.time.Day);
IS2datenum = convertTo(IS2dates,'yyyymmdd');
IS2dates_unique = unique(IS2dates);
IS2datenum_unique = convertTo(IS2dates_unique,'yyyymmdd');


SnowCoverfiles = dir([SnowCoverPath '/*.tif']);
SnowCoverdates = char(SnowCoverfiles.name);
SnowCoverdates = SnowCoverdates(:,[1:8]);
SnowCoverdates = str2num(SnowCoverdates);

%% calculate snow cover
for i = 1 %:length(IS2datenum)
    [d,ix] = min(abs(SnowCoverdates-IS2datenum_unique(i))); %find image date closest to IS2 date
    snowcover = imread([SnowCoverPath SnowCoverfiles(ix).name],'tif'); %read in the closest image
   
    IS2track = IS2(IS2datenum == IS2datenum_unique(i)); %index to only data on given date



end



























