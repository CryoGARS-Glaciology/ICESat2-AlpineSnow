function SnowCoverOut = SnowCover(IS2file,SnowCoverfilepath)
% Function to pull snow-on vs snow-free information from classified
% satilite imagry 
% 
% INPUTS: 
%       IS2file = csv of ICESat-2 data points
%       SnowCoverfilepath = path to folder with clasified snow cover images
% OUTPUTS:  
%       SnowCoverOut = 
% last modified Nov 2023 by Karina Zikan (karinazikan@u.boisestate.edu)

IS2dates = unique(datetime(IS2file.time.Year,IS2file.time.Month,IS2file.time.Day));
IS2datenum = convertTo(dates,'yyyymmdd');
IS2datenum = str2num(IS2datenum);

SnowCoverfiles = dir([SnowCoverPath '/*.tif']);
SnowCoverdates = char(SnowCoverfiles.name);
SnowCoverdates = SnowCoverdates(:,[1:8]);
SnowCoverdates = str2num(SnowCoverdates);

for i = 1:length(IS2datenum)
    [d,ix] = min(abs(SnowCoverdates-IS2datenum(i)));
    