%%% Make Snow Depth timeseries figure
%%%
%%% SPECIFIED INPUTS:


%% Inputs
clearvars;

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/snowex_toos_atl03_ai_20241121/';

%% Load data

%snotel data
granuel_files = dir([folderpath '*.csv']);
for j = 1:6
    granuel = readtable([folderpath  granuel_files(j+7).name]);
    % for i = 2:length(granuel_files)
    %     file = readtable([folderpath granuel_files(i).name]);
    %     granuel = cat(1,granuel,file);
    % end
    granuel(granuel.lon < -117,:) = []; granuel(granuel.lon > -115,:) = [];
    granuel(granuel.lat < 43,:) = []; granuel(granuel.lat > 44.5,:) = [];

    start_time = datetime(2018, 01, 01);
    time = start_time + seconds(granuel.delta_time);
    granuel.time = time; 

    %% group data
    G = findgroups(granuel.gt, granuel.yyyymmdd);
    group_ref_pa_deg = splitapply(@(x){x}, [granuel.ref_pa_deg] , G);
    group_time = splitapply(@(x){x}, [granuel.time] , G);

    %% plot
    % angle over time
    figure(j); clf
    hold on
    for i = 1:length(group_ref_pa_deg)
        plot(group_time{i}, group_ref_pa_deg{i}, 'LineWidth',3)
        xlabel('time'); ylabel('angle (deg)')
        set(gca,'fontsize',20);
    end
    hold off

    % location plot
    figure(j + 6); clf
    geoscatter(granuel.lat,granuel.lon,"filled")
    geobasemap topographic

end





