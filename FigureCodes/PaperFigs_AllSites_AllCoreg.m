%% Inputs
clearvars;
addpath(['/Users/karinazikan/Documents/ICESat2-AlpineSnow/functions'])
addpath(['/Users/karinazikan/Documents/cmocean'])
%colors
load('/Users/karinazikan/Documents/ScientificColourMaps8/vik/DiscretePalettes/vik10.mat');

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/';
%site abbreviation for file names
site_abbrevs = string({'RCEW';'Banner';'MCS';'DCEW'});
site_names = string({'RCEW';'BCS';'MCS';'DCEW'});
Coreg_type = string({'noCoreg';'Agg';'ByTrack'});

%Turn dtm or is2 slope correction
slope_correction = 2; % 0 = dtm, 1 = is2, 2 = no slope correction

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
    for k = 1:length(Coreg_type)
        filepath = strcat(folderpath, site_abbrevs(j), '/IS2_Data/A6-40/ATL06-A6-40-AllData-', Coreg_type(k), '_on.csv');
        df_on{j,k} = readtable(filepath);
        filepath = strcat(folderpath, site_abbrevs(j), '/IS2_Data/A6-40/ATL06-A6-40-AllData-', Coreg_type(k), '_off.csv');
        df_off{j,k} = readtable(filepath);
        % slope correction
        if slope_correction == 1
            df_off{j,k}.elev_residuals_vertcoreg = df_off{j,k}.elev_residuals_vertcoreg_is2_slopecorrected;
            df_on{j,k}.elev_residuals_vertcoreg = df_on{j,k}.elev_residuals_vertcoreg_is2_slopecorrected;
            df_off{j,k}.slope_mean = df_off{j,k}.IS2_slope_deg;
            df_on{j,k}.slope_mean = df_on{j,k}.IS2_slope_deg;
            if k ==1 && j == 1
                ('IS2 Slope correction used')
            else
            end
        elseif slope_correction == 2
            if k ==1 && j == 1
                ('No Slope correction used')
            else
            end
        else
            df_off{j,k}.elev_residuals_vertcoreg = df_off{j,k}.elev_residuals_vertcoreg_dtm_slopecorrected;
            df_on{j,k}.elev_residuals_vertcoreg = df_on{j,k}.elev_residuals_vertcoreg_dtm_slopecorrected;
            df_off{j,k}.slope_mean = df_off{j,k}.IS2_slope_deg;
            df_on{j,k}.slope_mean = df_on{j,k}.IS2_slope_deg;
            if k ==1 && j == 1
                ('DTM Slope correction used')
            else
            end
        end
    end

    %snotel data
    snotel_files = dir(strcat(folderpath, site_abbrevs(j), '/snotel/*.csv'));
    snotel{j} = readtable(strcat(folderpath, site_abbrevs(j), '/snotel/', snotel_files(1).name));
    for i = 2:length(snotel_files)
        file = readtable(strcat(folderpath, site_abbrevs(j), '/snotel/', snotel_files(i).name));
        snotel{j} = cat(1,snotel{j},file);
    end
    snotel{j}.SNWD_I_1_in_ = snotel{j}.SNWD_I_1_in_ * 0.0254; %convert from in to m
    snotel{j}.SNWD_I_1_in_(snotel{j}.SNWD_I_1_in_ < 0) = NaN;
end


%% Plots
%% Snow depth time series
%timeseries
for k = 1:length(Coreg_type)
    % figure(k); clf
    % for j = 1:length(site_abbrevs)
    %     df_on_temp = df_on{j,k};
    %     df_off_temp = df_off{j,k};
    %     snotel_temp = snotel{j};
    %
    %     snowdepth = table([datetime(df_on_temp.time.Year,df_on_temp.time.Month,df_on_temp.time.Day)], df_on_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
    %     snowdepth_dategroup = varfun(@(x)median(x,'omitnan'),snowdepth,'GroupingVariables',{'time'});
    %     elevresiduals = table([datetime(df_off_temp.time.Year,df_off_temp.time.Month,df_off_temp.time.Day)], df_off_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
    %     elevresiduals_dategroup = varfun(@(x)median(x,'omitnan'),elevresiduals,'GroupingVariables',{'time'});
    %     snoteldepth = table(snotel_temp.Date, snotel_temp.SNWD_I_1_in_ ,'VariableNames',["time","SnowDepth"]);
    %     snoteldepth_dategroup = varfun(@(x)nanmean(x),snoteldepth,'GroupingVariables',{'time'});
    %
    %     subplot(4,1,j); hold on
    %     plot(snoteldepth_dategroup.time,snoteldepth_dategroup.Fun_SnowDepth, 'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]	)
    %     scatter(snowdepth_dategroup.time,snowdepth_dategroup.Fun_residuals,50,'filled','MarkerFaceColor',[0, 0.4470, 0.7410])
    %     title(site_names{j});
    %     set(gca,'fontsize',20);
    %     ylabel('Snow Depth (m)')
    %     hold off
    % end
    % legend('SNOTEL mean snow depth','ICESat-2 median snow depth')
    % sgtitle(Coreg_type{k});
    %
    % figure(k+3); clf
    % for j = 1:length(site_abbrevs)
    %     df_on_temp = df_on{j,k};
    %     df_off_temp = df_off{j,k};
    %     snotel_temp = snotel{j};
    %
    %     snowdepth = table([datetime(df_on_temp.time.Year,df_on_temp.time.Month,df_on_temp.time.Day)], df_on_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
    %     snowdepth_dategroup = varfun(@(x)median(x,'omitnan'),snowdepth,'GroupingVariables',{'time'});
    %     elevresiduals = table([datetime(df_off_temp.time.Year,df_off_temp.time.Month,df_off_temp.time.Day)], df_off_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
    %     elevresiduals_dategroup = varfun(@(x)median(x,'omitnan'),elevresiduals,'GroupingVariables',{'time'});
    %     snoteldepth = table(snotel_temp.Date, snotel_temp.SNWD_I_1_in_ ,'VariableNames',["time","SnowDepth"]);
    %     snoteldepth_dategroup = varfun(@(x)nanmean(x),snoteldepth,'GroupingVariables',{'time'});
    %
    %     subplot(4,1,j); hold on
    %     plot(snoteldepth_dategroup.time,snoteldepth_dategroup.Fun_SnowDepth, 'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]	)
    %     scatter(elevresiduals_dategroup.time,elevresiduals_dategroup.Fun_residuals,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
    %     scatter(snowdepth_dategroup.time,snowdepth_dategroup.Fun_residuals,'filled','MarkerFaceColor',[0, 0.4470, 0.7410])
    %     set(gca,'fontsize',20);
    %     title(site_names{j});
    %     ylabel('Snow Depth (m)')
    %     hold off
    % end
    % legend('SNOTEL mean snow depth','ICESat-2 median snow free elevation residuals','ICESat-2 median snow depth')
    % sgtitle(Coreg_type{k});

    figure(k+6); clf
    for j = 1:length(site_abbrevs)
        df_on_temp = df_on{j,k};
        df_off_temp = df_off{j,k};
        snotel_temp = snotel{j};

        if j == 1
            df_off_temp(df_off_temp.time.Month > 11,:) = [];
            df_off_temp(df_off_temp.time.Month < 5,:) = [];
        else
            df_off_temp(df_off_temp.time.Month > 10,:) = [];
            df_off_temp(df_off_temp.time.Month < 6,:) = [];
        end

        snowdepth = table([datetime(df_on_temp.time.Year,df_on_temp.time.Month,df_on_temp.time.Day)], df_on_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        snowdepth_dategroup = varfun(@(x)median(x,'omitnan'),snowdepth,'GroupingVariables',{'time'});
        elevresiduals = table([datetime(df_off_temp.time.Year,df_off_temp.time.Month,df_off_temp.time.Day)], df_off_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        elevresiduals_dategroup = varfun(@(x)median(x,'omitnan'),elevresiduals,'GroupingVariables',{'time'});
        elevresiduals_error_dategroup = varfun(@(x)iqr(x),elevresiduals,'GroupingVariables',{'time'});
        snoteldepth = table(snotel_temp.Date, snotel_temp.SNWD_I_1_in_ ,'VariableNames',["time","SnowDepth"]);
        snoteldepth_dategroup = varfun(@(x)nanmean(x),snoteldepth,'GroupingVariables',{'time'});

        subplot(4,1,j); cla; hold on
        %plot(snoteldepth_dategroup.time,snoteldepth_dategroup.Fun_SnowDepth, 'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]	)
        %errorbar(elevresiduals_dategroup.time,elevresiduals_dategroup.Fun_residuals,elevresiduals_error_dategroup.Fun_residuals, 'vertical', 'LineStyle', 'none');
        scatter(elevresiduals_dategroup.time,elevresiduals_dategroup.Fun_residuals,'filled','o')
        %scatter(snowdepth_dategroup.time,snowdepth_dategroup.Fun_residuals,'filled','MarkerFaceColor',[0, 0.4470, 0.7410])
        yline(0)
        set(gca,'fontsize',20);
        title(site_names{j});
        ylabel('Elevation Residual (m)')
        xlim([datetime(2018,06,01) datetime(2024,01,01)])
        ylim([-1 1])
        hold off
    end
    legend( 'ICESat-2 median snow free elevation residuals','zero line')
    sgtitle(Coreg_type{k});
end

