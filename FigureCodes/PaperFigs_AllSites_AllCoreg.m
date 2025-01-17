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
slope_correction = 1; % 0 = dtm, 1 = is2, 2 = no slope correction
%Turn dtm or is2 slope correction
slope_filter = 1; % 0 = none, 1 = remove slopes > 30 degrees
% ICESat-2 residuals below 0 to NaN?
remove_negative = 0; % 0 = off, 1 = on

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
        end
        if k ==1 && j == 1
            ('DTM Slope correction used')
        else
        end
        if slope_filter == 1
            df_off{j,k}(df_off{j,k}.slope_mean > 30,:) = [];
            df_on{j,k}(df_on{j,k}.slope_mean > 30,:) = [];
        else
        end
        % remove negative ICESat-2 elevation residuals
        if remove_negative == 1
            df_on{j,k}.elev_residuals_vertcoreg(df_on{j,k}.elev_residuals_vertcoreg < 0) = NaN;
            if k ==1 && j == 1
                ('Negative ICESat-2 snow depth removed')
            else
            end
        else
        end
    end


    %snotel data
    snotel_files = dir(strcat(folderpath, site_abbrevs(j), '/snotel/*.csv'));
    if j == 1 %RCEW site
        temp_file = readtable(strcat(folderpath, site_abbrevs(j), '/snotel/', snotel_files(1).name));
        Dates = datetime(temp_file.date_time.Year,temp_file.date_time.Month,temp_file.date_time.Day);
        snotel{j}.Date = [Dates; Dates; Dates; Dates];
        snotel{j}.SNWD_I_1_in_ = [temp_file.RME_176; temp_file.RME_176b; temp_file.RME_rmsp3; temp_file.RME_rmsp3b];
        snotel{j}.SNWD_I_1_in_ = snotel{j}.SNWD_I_1_in_ / 100;
    else
        snotel{j} = readtable(strcat(folderpath, site_abbrevs(j), '/snotel/', snotel_files(1).name));
        for i = 2:length(snotel_files)
            file = readtable(strcat(folderpath, site_abbrevs(j), '/snotel/', snotel_files(i).name));
            snotel{j} = cat(1,snotel{j},file);
        end
        snotel{j}.SNWD_I_1_in_ = snotel{j}.SNWD_I_1_in_ * 0.0254; %convert from in to m
        snotel{j}.SNWD_I_1_in_(snotel{j}.SNWD_I_1_in_ < 0) = NaN;
    end
end

%% Stats
for k = 1:length(Coreg_type)
    for j = 1:length(site_abbrevs)
        Rnmad(j,k) = 1.4826*median(abs(df_off{j,k}.elev_residuals_vertcoreg-nanmean(df_off{j,k}.elev_residuals_vertcoreg)),'omitnan'); % normalized meadian absolute difference
        RMSE(j,k) = sqrt(nanmean((df_off{j,k}.elev_residuals_vertcoreg).^2)); % Root mean square error
    end
end

%% Plots
%% Coregistration types

    figure(12); clf
    xbounds = [-1,1];
    for  j = 1:length(site_abbrevs)
        if j < 3
            xbounds = [-2,2];
        else
            xbounds = [-3,3];
        end
        subplot(2,2,j); hold on
        for k = 1:length(Coreg_type)
            pd = fitdist(df_off{j,k}.elev_residuals_vertcoreg,'kernel','Kernel','normal');
            fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 3);
        end
        title(site_names{j});
        xline(0, 'Linewidth', 1,'Color','black');
        set(gca,'fontsize',20,'xlim',xbounds);
        xlabel('Snow free Vertical offset (m)'); ylabel('Probability density');
    end
    %set(gcf,'position',[50 50 800 400]);
    
    legend('No Coregistration','Aggregated Coregistration','By Track Coregistration');
    hold off


%% Snow depth time series
%timeseries
for k = 1:length(Coreg_type)
    fig = figure(k); clf;
    for j = 1:length(site_abbrevs)
        df_on_temp = df_on{j,k};
        df_off_temp = df_off{j,k};
        snotel_temp = snotel{j};

        snowdepth = table([datetime(df_on_temp.time.Year,df_on_temp.time.Month,df_on_temp.time.Day)], df_on_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        snowdepth_dategroup = varfun(@(x)median(x,'omitnan'),snowdepth,'GroupingVariables',{'time'});
        elevresiduals = table([datetime(df_off_temp.time.Year,df_off_temp.time.Month,df_off_temp.time.Day)], df_off_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        elevresiduals_dategroup = varfun(@(x)median(x,'omitnan'),elevresiduals,'GroupingVariables',{'time'});
        snoteldepth = table(snotel_temp.Date, snotel_temp.SNWD_I_1_in_ ,'VariableNames',["time","SnowDepth"]);
        snoteldepth_dategroup = varfun(@(x)nanmean(x),snoteldepth,'GroupingVariables',{'time'});
        snoteldepth_dategroup_IS2_dates = snoteldepth_dategroup(ismember(snoteldepth_dategroup.time, elevresiduals_dategroup.time), :);

        subplot(4,1,j); hold on
        plot(snoteldepth_dategroup.time,snoteldepth_dategroup.Fun_SnowDepth, 'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250]	)
        scatter(snowdepth_dategroup.time,snowdepth_dategroup.Fun_residuals,75,'filled','MarkerFaceColor',[0, 0.4470, 0.7410])
        title(site_names{j});
        set(gca,'fontsize',20);
        xlim([datetime(2018,1,1) datetime(2024,6,1)])
        ylim([-0.5 4])
        ylabel('Snow Depth')
        
        hold off
    end
    xlabel('Date'); ylabel('Snow Depth')
    legend('SNOTEL mean snow depth','ICESat-2 median snow depth')
    sgtitle(Coreg_type{k});
    % han=axes(fig,'visible','off');
    % han.Title.Visible='on';
    % han.XLabel.Visible='on';
    % han.YLabel.Visible='on';
    % set(han,'fontsize',20);
    % ylabel(han,'Snow Depth (m)');
    % xlabel(han,'Date');
    % title(han,Coreg_type{k});
    % han.TitleHorizontalAlignment = 'left';

    figure(k+3); clf
    for j = 1:length(site_abbrevs)
        df_on_temp = df_on{j,k};
        df_off_temp = df_off{j,k};
        snotel_temp = snotel{j};

        snowdepth = table([datetime(df_on_temp.time.Year,df_on_temp.time.Month,df_on_temp.time.Day)], df_on_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        snowdepth_dategroup = varfun(@(x)median(x,'omitnan'),snowdepth,'GroupingVariables',{'time'});
        elevresiduals = table([datetime(df_off_temp.time.Year,df_off_temp.time.Month,df_off_temp.time.Day)], df_off_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        elevresiduals_dategroup = varfun(@(x)median(x,'omitnan'),elevresiduals,'GroupingVariables',{'time'});
        snoteldepth = table(snotel_temp.Date, snotel_temp.SNWD_I_1_in_ ,'VariableNames',["time","SnowDepth"]);
        snoteldepth_dategroup = varfun(@(x)nanmean(x),snoteldepth,'GroupingVariables',{'time'});

        subplot(4,1,j); hold on
        plot(snoteldepth_dategroup.time,snoteldepth_dategroup.Fun_SnowDepth, 'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]	)
        scatter(elevresiduals_dategroup.time,elevresiduals_dategroup.Fun_residuals,75,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
        scatter(snowdepth_dategroup.time,snowdepth_dategroup.Fun_residuals,75,'filled','MarkerFaceColor',[0, 0.4470, 0.7410])
        set(gca,'fontsize',20);
        title(site_names{j});
        hold off
        xlim([datetime(2018,1,1) datetime(2024,6,1)])
        ylim([-0.5 4])
        ylabel('Snow Depth')
    end
    legend('AWS mean snow depth','ICESat-2 median snow free elevation residuals','ICESat-2 median snow depth')
    sgtitle(Coreg_type{k});
     xlabel('Date');
    % han=axes(fig,'visible','off');
    % han.Title.Visible='on';
    % han.XLabel.Visible='on';
    % han.YLabel.Visible='on';
    % set(han,'fontsize',20);
    % ylabel(han,'Snow Depth (m)');
    % xlabel(han,'Date');
    % title(han,Coreg_type{k});
    % han.TitleHorizontalAlignment = 'left';

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
        scatter(elevresiduals_dategroup.time,elevresiduals_dategroup.Fun_residuals,75,'filled','o')
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
%%
fig = figure(10); clf
for j = 1:length(site_abbrevs)
    snotel_temp = snotel{j};
    snoteldepth = table(snotel_temp.Date, snotel_temp.SNWD_I_1_in_ ,'VariableNames',["time","SnowDepth"]);
    snoteldepth_dategroup = varfun(@(x)nanmean(x),snoteldepth,'GroupingVariables',{'time'});

    subplot(4,1,j); hold on
    plot(snoteldepth_dategroup.time,snoteldepth_dategroup.Fun_SnowDepth, 'LineWidth',2)

    for k = 1:length(Coreg_type)
        df_on_temp = df_on{j,k};
        df_off_temp = df_off{j,k};
        snotel_temp = snotel{j};

        snowdepth = table([datetime(df_on_temp.time.Year,df_on_temp.time.Month,df_on_temp.time.Day)], df_on_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        snowdepth_dategroup = varfun(@(x)median(x,'omitnan'),snowdepth,'GroupingVariables',{'time'});
        elevresiduals = table([datetime(df_off_temp.time.Year,df_off_temp.time.Month,df_off_temp.time.Day)], df_off_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        elevresiduals_dategroup = varfun(@(x)median(x,'omitnan'),elevresiduals,'GroupingVariables',{'time'});

        scatter(snowdepth_dategroup.time,snowdepth_dategroup.Fun_residuals,75,'filled')
        title(site_names{j});
        set(gca,'fontsize',20);
    end
    hold off
end
legend('AWS mean snow depth','no coreg', 'agg', 'bytrack')
%%
for k = 2
figure(11); clf;
    for j = 1:4%length(site_abbrevs)
        df_on_temp = df_on{j,k};
        df_off_temp = df_off{j,k};
        snotel_temp = snotel{j};

        snowdepth = table([datetime(df_on_temp.time.Year,df_on_temp.time.Month,df_on_temp.time.Day)], df_on_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        snowdepth_dategroup = varfun(@(x)median(x,'omitnan'),snowdepth,'GroupingVariables',{'time'});
        elevresiduals = table([datetime(df_off_temp.time.Year,df_off_temp.time.Month,df_off_temp.time.Day)], df_off_temp.elev_residuals_vertcoreg, 'VariableNames',["time","residuals"]);
        elevresiduals_dategroup = varfun(@(x)median(x,'omitnan'),elevresiduals,'GroupingVariables',{'time'});
        snoteldepth = table(snotel_temp.Date, snotel_temp.SNWD_I_1_in_ ,'VariableNames',["time","SnowDepth"]);
        snoteldepth_dategroup = varfun(@(x)nanmean(x),snoteldepth,'GroupingVariables',{'time'});
        snoteldepth_dategroup_IS2_dates = snoteldepth_dategroup(ismember(snoteldepth_dategroup.time, elevresiduals_dategroup.time), :);

        subplot(4,1,j); hold on
        scatter(elevresiduals_dategroup.GroupCount(1:height(snoteldepth_dategroup_IS2_dates)),(elevresiduals_dategroup.Fun_residuals(1:height(snoteldepth_dategroup_IS2_dates))-snoteldepth_dategroup_IS2_dates.Fun_SnowDepth),75,'filled','MarkerFaceColor',[0, 0.4470, 0.7410])
        yline(0)
        title(site_names{j});
        set(gca,'fontsize',20);
        if j == 2
            ylabel('ICESat-2 Snow Depth - AWS Snow Depth')
        end
     
        
        hold off
    end
    xlabel('# of ICESat-2 observations'); 

end
% han=axes(fig,'visible','off');
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% set(han,'fontsize',20);
% ylabel(han,'Snow Depth (m)');
% xlabel(han,'Date');
%title(han,Coreg_type{k});
%   han.TitleHorizontalAlignment = 'left';
