%% Inputs
clearvars;
addpath('/Users/karinazikan/Documents/ICESat2-AlpineSnow/functions')
addpath('/Users/karinazikan/Documents/cmocean')

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/';
%site abbreviation for file names
abbrev = 'MCS';

outputname = [ abbrev '_segments_and_gtx_date_first'];

%% Load data ({1} = ATL08, {2} = ATL06, {3} = ATL06 w/ ATL08 classification)
%load the reference elevation data
filepath = strcat(folderpath, abbrev, '/IS2_Data/A6-40/ATL06-A6-40-AllData-Agg');
df = readtable(filepath);
df.time = datetime(df.time.Year,df.time.Month,df.time.Day);
df.Easting = df.Easting./(10^3);
df.Northing = df.Northing./(10^3);
%snow on + snow off
df_off = df(df.snowcover == 0, :);
df_on = df(df.snowcover == 1, :);

%unique dates
unique_dates = unique(df.time);
% 
% %group segment ids and gtx
% groups = groupsummary(df,["segment_id","gt"]);
% %remame gts
% gt(groups.gt == 10,:) = "gt1l";
% gt(groups.gt == 20,:) = "gt1r";
% gt(groups.gt == 30,:) = "gt2l";
% gt(groups.gt == 40,:) = "gt2r";
% gt(groups.gt == 50,:) = "gt3l";
% gt(groups.gt == 60,:) = "gt3r";
% 
% Segmesnts_gtx = [groups.segment_id, gt];
Segmesnts_gtx = [];

for i = 1:length(unique_dates)
    df_loop = df(df.time == unique_dates(i),:);
    %remame gts
    gt(df_loop.gt == 10,:) = "gt1l";
    gt(df_loop.gt == 20,:) = "gt1r";
    gt(df_loop.gt == 30,:) = "gt2l";
    gt(df_loop.gt == 40,:) = "gt2r";
    gt(df_loop.gt == 50,:) = "gt3l";
    gt(df_loop.gt == 60,:) = "gt3r";
    Segmesnts_gtx = [Segmesnts_gtx; df_loop.segment_id(1), gt(1)];
    clear gt
end

%% export
writematrix(Segmesnts_gtx,[folderpath abbrev '/IS2_Data/' outputname '.csv']);

