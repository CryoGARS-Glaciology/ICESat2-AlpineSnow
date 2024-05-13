%% Inputs
clearvars;
addpath(['/Users/karinazikan/Documents/ICESat2-AlpineSnow/functions'])
addpath(['/Users/karinazikan/Documents/cmocean'])

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/';
%site abbreviation for file names
site_abbrevs = string({'RCEW';'Banner';'MCS';'DCEW'});
site_names = string({'RCEW';'BCS';'MCS';'DCEW'});

%Turn dtm or is2 slope correction
slope_correction = 0; % 0 = dtm, 1 = is2, 2 = no correction

slope_cutoff = 20;

%Weather Station Location
% Reynolds snowex: snotel_E = 519729; snotel_N = 4768225;
% Banner snotel: snotel_E = 640823; snotel_N = 4907084;
% MCS snotel: snotel_E = 607075; snotel_N = 4865193;
% DCEW little deer point AWS: snotel_E = 570697; snotel_N = 4843042;
snotel_E = {519729; 640823; 607075; 570697};
snotel_N = {4768225; 4907084; 4865193; 4843042};

%set colors
colors{1} = cmocean('-dense',6);
colors{2} = cmocean('-algae',5);
colors{3} = cmocean('ice',5);
colors{4} = cmocean('-amp',5);

%% Read in data
for j = 1:length(site_abbrevs)
    % IS2 data
    filepath = strcat(folderpath, site_abbrevs(j), '/IS2_Data/ATL06-atl08class-AllData_on.csv');
    df_on{j} = readtable(filepath);
    filepath = strcat(folderpath, site_abbrevs(j), '/IS2_Data/ATL06-atl08class-AllData_off.csv');
    df_off{j} = readtable(filepath);

    %snotel data
    snotel_files = dir(strcat(folderpath, site_abbrevs(j), '/snotel/*.csv'));
    snotel{j} = readtable(strcat(folderpath, site_abbrevs(j), '/snotel/', snotel_files(1).name));
    for i = 2:length(snotel_files)
        file = readtable(strcat(folderpath, site_abbrevs(j), '/snotel/', snotel_files(i).name));
        snotel{j} = cat(1,snotel{j},file);
    end
    snotel{j}.SNWD_I_1_in_ = snotel{j}.SNWD_I_1_in_ * 0.0254; %convert from in to m
    snotel{j}.SNWD_I_1_in_(snotel{j}.SNWD_I_1_in_ < 0) = NaN;

    %%
    if slope_correction == 1
        df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected = df_off{j}.elev_residuals_vertcoreg_is2_slopecorrected;
        df_on{j}.elev_residuals_vertcoreg_dtm_slopecorrected = df_on{j}.elev_residuals_vertcoreg_is2_slopecorrected;
        df_off{j}.slope_mean = df_off{j}.IS2_slope_deg;
        df_on{j}.slope_mean = df_on{j}.IS2_slope_deg;
        ('IS2 Slope correction used')
    elseif slope_correction == 2
        df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected = df_off{j}.elev_residuals_vertcoreg;
        df_on{j}.elev_residuals_vertcoreg_dtm_slopecorrected = df_on{j}.elev_residuals_vertcoreg;
        ('No Slope correction used')
    else
        ('DTM Slope correction used')
    end
end

%% Site stats
for j = 1:length(site_abbrevs)
    med(j) = median(df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected,'omitnan');
    nmad(j) = 1.4826*median(abs(df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected-nanmean(df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected)),'omitnan');
    %below slope cutoff
    ix = find(df_off{j}.IS2_slope_deg < slope_cutoff);
    median_below(j) = median(df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected(ix),'omitnan');
    nmad_below(j) = 1.4826*median(abs(df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected(ix)-nanmean(df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected(ix))),'omitnan');
    %above slope cutoff
    ix = find(df_off{j}.IS2_slope_deg > slope_cutoff);
    median_above(j) = median(df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected(ix),'omitnan');
    nmad_above(j) = 1.4826*median(abs(df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected(ix)-nanmean(df_off{j}.elev_residuals_vertcoreg_dtm_slopecorrected(ix))),'omitnan');   
end