%% Inputs
clearvars;
addpath(['/Users/karinazikan/Documents/ICESat2-AlpineSnow/functions'])
addpath(['/Users/karinazikan/Documents/cmocean'])

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/';
%site abbreviation for file names
site_abbrevs = string({'RCEW';'Banner';'MCS';'DCEW_2'});
site_names = string({'RCEW';'BCS';'MCS';'DCEW'});

%Turn slope correction off or on
slope_correction = 0; % 0 = off, 1 = on

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
        %calculate quadratic slope correction
        slope = df_off{j}.slope_mean;
        x= slope; y = df_off{j}.elev_residuals_vertcoreg;
        ind = isnan(x) | isnan(y); %index nans
        x(ind) = []; y(ind) = []; %remove nans
        p = polyfit(x,y,2); % fit quadratic
        % Vertical corregistration
        df_on{j}.elev_residuals_vertcoreg = df_on{j}.elev_residuals_vertcoreg-polyval(p,df_on{j}.slope_mean);
        df_off{j}.elev_residuals_vertcoreg = df_off{j}.elev_residuals_vertcoreg-polyval(p,df_off{j}.slope_mean);
        ('Slope correction applied')
    else
        ('Slope correction not applied')
    end
end


%
%
%
%% Figures
%% Non-parametric pdfs
% snow on vs snow pff
figure(1); clf; hold on
for j = 1:length(site_abbrevs)
    subplot(length(site_abbrevs),1,j); hold on
    pd = fitdist(df_on{j}.elev_residuals_vertcoreg,'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2);
    pd = fitdist(df_off{j}.elev_residuals_vertcoreg,'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2);
    set(gca,'fontsize',16,'xlim',[-3 3]);
    %set(gcf,'position',[50 50 800 400]);
    xlabel('Snow free Vertical offset (m)'); ylabel('Probability density');
    legend('Snow On','Snow Off');
    title(site_names{j});
end
hold off

figure(2); clf; hold on
for j = 1:length(site_abbrevs)
    subplot(length(site_abbrevs),1,j); hold on
    pd = fitdist(df_on{j}.elev_residuals_vertcoreg,'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2);
    pd = fitdist(df_off{j}.elev_residuals_vertcoreg,'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2);
    pd = fitdist(df_on{j}.elev_residuals_vertcoreg_slopecorrected,'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2);
    pd = fitdist(df_off{j}.elev_residuals_vertcoreg_slopecorrected,'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2);
    set(gca,'fontsize',16,'xlim',[-3 3]);
    %set(gcf,'position',[50 50 800 400]);
    xlabel('Snow free Vertical offset (m)'); ylabel('Probability density');
    legend('Snow On','Snow Off','Slope Corrected Snow On','Slope Corrected Snow Off');
    title(site_names{j});
end
hold off

%% Terrain Histograms & boxplots
Nbin = 8; %set humber of bins for box plots
figure(3); clf
for j = 1:length(site_abbrevs)
    %ELEVATION
    subplot(length(site_abbrevs),3,(3*j-2));
    h = histogram([df_off{j}.elevation_report_mean; df_on{j}.elevation_report_mean],Nbin,'FaceColor',colors{1}(3,:));
    elev_binwidth{j} = h.BinWidth; elev_binedges{j} = h.BinEdges;
    clear h;
    xlabel('Elevation (m a.s.l.)','fontsize',16); ylabel('Observations','fontsize',16);
    %ASPECT
    subplot(length(site_abbrevs),3,(3*j-1));
    h = histogram([df_off{j}.aspect_mean; df_on{j}.aspect_mean],Nbin,'FaceColor',colors{2}(3,:));
    aspect_binwidth{j} = h.BinWidth; aspect_binedges{j} = h.BinEdges;
    clear h;
    xlabel('Aspect (degrees)','fontsize',16); 
    title(site_names{j},'fontsize',16);
    %SLOPE
    subplot(length(site_abbrevs),3,(3*j));
    h = histogram([df_off{j}.slope_mean; df_on{j}.slope_mean],Nbin,'FaceColor',colors{3}(3,:));
    slope_binwidth{j} = h.BinWidth; slope_binedges{j} = h.BinEdges;
    clear h;
    xlabel('Slope (degrees)','fontsize',16);
end

% Terrain Boxplots
% no slope correction
figure(4); clf
for j = 1:length(site_abbrevs)
        %ELEVATION
    bins = {num2str(elev_binedges{j}(2))};
    for i= 3:length(elev_binedges{j})
        bins = [bins; {num2str(elev_binedges{j}(i))}];
    end
    subplot(length(site_abbrevs),3,(3*j-2));
    hold on
    groupElev = discretize(df_off{j}.elevation_report_mean,elev_binedges{j},'categorical',bins);
    boxchart(groupElev, df_off{j}.elev_residuals_vertcoreg,'BoxFaceColor',colors{1}(3,:),'MarkerStyle','none')
    ylim([-2,2])
    set(gca,'fontsize',16,'box','on'); drawnow;
    xlabel('Elevation (m a.s.l.)','fontsize',16); ylabel('Elevation residuals (m)','fontsize',16);
   
        %ASPECT
    bins = {num2str(aspect_binedges{j}(2))};
    for i= 3:length(aspect_binedges{j})
        bins = [bins; {num2str(aspect_binedges{j}(i))}];
    end
    subplot(length(site_abbrevs),3,(3*j-1));
    hold on
    groupAspect = discretize(df_off{j}.aspect_mean,aspect_binedges{j},'categorical',bins);
    boxchart(groupAspect, df_off{j}.elev_residuals_vertcoreg,'BoxFaceColor',colors{2}(3,:),'MarkerStyle','none')
    ylim([-2,2])
    set(gca,'fontsize',16,'box','on'); drawnow;
    xlabel('Aspect (degrees)','fontsize',16); 
    title(site_names{j},'fontsize',16);
    
        %SLOPE
    bins = {num2str(slope_binedges{j}(2))};
    for i= 3:length(slope_binedges{j})
        bins = [bins; {num2str(slope_binedges{j}(i))}];
    end
    subplot(length(site_abbrevs),3,(3*j));
    hold on
    groupSlope = discretize(df_off{j}.slope_mean,slope_binedges{j},'categorical',bins);
    boxchart(groupSlope, df_off{j}.elev_residuals_vertcoreg,'BoxFaceColor',colors{3}(3,:),'MarkerStyle','none')
    ylim([-3,3])
    set(gca,'fontsize',16,'box','on'); drawnow;
    xlabel('Slope (degrees)','fontsize',16); 
end

% slope correction
figure(5); clf
for j = 1:length(site_abbrevs)
    %ELEVATION
    bins = {num2str(elev_binedges{j}(2))};
    for i= 3:length(elev_binedges{j})
        bins = [bins; {num2str(elev_binedges{j}(i))}];
    end
    subplot(length(site_abbrevs),3,(3*j-2));
    hold on
    groupElev = discretize(df_off{j}.elevation_report_mean,elev_binedges{j},'categorical',bins);
    boxchart(groupElev, df_off{j}.elev_residuals_vertcoreg_slopecorrected,'BoxFaceColor',colors{1}(3,:),'MarkerStyle','none')
    ylim([-2,2])
    set(gca,'fontsize',16,'box','on'); drawnow;
    xlabel('Elevation (m a.s.l.)','fontsize',16); ylabel('Elevation residuals (m)','fontsize',16);
   
    %ASPECT
    bins = {num2str(aspect_binedges{j}(2))};
    for i= 3:length(aspect_binedges{j})
        bins = [bins; {num2str(aspect_binedges{j}(i))}];
    end
    subplot(length(site_abbrevs),3,(3*j-1));
    hold on
    groupAspect = discretize(df_off{j}.aspect_mean,aspect_binedges{j},'categorical',bins);
    boxchart(groupAspect, df_off{j}.elev_residuals_vertcoreg_slopecorrected,'BoxFaceColor',colors{2}(3,:),'MarkerStyle','none')
    ylim([-2,2])
    set(gca,'fontsize',16,'box','on'); drawnow;
    xlabel('Aspect (degrees)','fontsize',16); 
    title(site_names{j},'fontsize',16);
    
    %SLOPE
    bins = {num2str(slope_binedges{j}(2))};
    for i= 3:length(slope_binedges{j})
        bins = [bins; {num2str(slope_binedges{j}(i))}];
    end
    subplot(length(site_abbrevs),3,(3*j));
    hold on
    groupSlope = discretize(df_off{j}.slope_mean,slope_binedges{j},'categorical',bins);
    boxchart(groupSlope, df_off{j}.elev_residuals_vertcoreg_slopecorrected,'BoxFaceColor',colors{3}(3,:),'MarkerStyle','none')
    ylim([-3,3])
    set(gca,'fontsize',16,'box','on'); drawnow;
    xlabel('Slope (degrees)','fontsize',16); 
end
