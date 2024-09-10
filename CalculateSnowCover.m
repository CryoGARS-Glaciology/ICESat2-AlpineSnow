%%% This code was writen by Karina Zikan and Ellyn Enderlin to pull snow-on 
%%% snow-free information from classified NDSI images
%%%
%%% SPECIFIED INPUTS:
%%%     folderpath = path to site folder
%%%     IS2file = icesat2 csv filename
%%%     abbrev = site abreviation
%%%     acronym = ICESat-2 product acronym (ATL08, ALT06, ect)
%%%     outputname = output file name
%%%
%%% OUTPUTS:
%%%     IS2 = icesat2 csv with snow cove and reference snow cover map dates
%%%     appened
%%%
%%% Last updated: Nov 2023 by Karina Zikan

%% Inputs
clearvars; close all;

%Folder path 
folderpath = '/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/RCEW/';
% folderpath = '/Users/ellynenderlin/Research/NASA_CryoIdaho/mountains/RCEW/snow-mapping/';
cd(folderpath);
SnowCoverPath = [folderpath 'SnowCover/'];
%site abbreviation for file names
abbrev = 'RCEW';
% 
% % Site shapefile 
shapefile = 'RCEW-outline_WGS84.shp';

%filename
IS2file = 'RCEW-ICESat2-ATL06-atl08class.csv';
acronym = 'ATL06';

%Set output name
outputname = [abbrev,'-ICESat2-',acronym,'-atl08class-SnowCover.csv'];

%% read in ICESat-2 data
icesat2path = [folderpath 'IS2_Data/' IS2file];
IS2 = readtable(icesat2path);

watershed_border = shaperead([folderpath 'ROIs/' shapefile]);

%% Set up date arrays
IS2dates = datetime(IS2.time.Year,IS2.time.Month,IS2.time.Day);
IS2datenum = convertTo(IS2dates,'yyyymmdd');
IS2dates_unique = unique(IS2dates);
IS2datenum_unique = convertTo(IS2dates_unique,'yyyymmdd');


SnowCoverfiles = dir([SnowCoverPath '/*.tif']);
% SnowCoverdates = char(SnowCoverfiles.name);
% SnowCoverdates = SnowCoverdates(:,[1:8]);
SnowCoverdates = char(SnowCoverfiles.name);
SnowCoverdates = SnowCoverdates(:,[1:8]);
SnowCoverdates = str2num(SnowCoverdates);

%% calculate snow cover
snowvector = NaN(size(IS2datenum)); %set up a dummy vector to hold the binary snow flag to add to the table
snowdatevector = zeros(size(IS2datenum));
for i = 1:size(IS2datenum_unique,1) %8
    disp(num2str(IS2datenum_unique(i)));
    [d,ix] = min(abs(SnowCoverdates-IS2datenum_unique(i))); %find image date closest to IS2 date
%     snowcover = imread([SnowCoverPath SnowCoverfiles(ix).name],'tif'); %read in the closest image
    [snowcover,R] = readgeoraster([SnowCoverPath SnowCoverfiles(ix).name]); %read in the closest image
    snowX = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
    snowY = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
    [snowXgrid,snowYgrid] = meshgrid(snowX,snowY); clear R;
    
    disp(SnowCoverfiles(ix).name);

    %ONLY NEED FOR TESTING WITH IMAGES, NOT SNOW MAPS, ASSUMING snowcover =
    %1 if there is snow and snowcover = 0 if there is no snow
    % % disp(SnowCoverfiles(ix).name);
    % % if contains(SnowCoverfiles(ix).name,'Sentinel')
    % %     NDSI = (snowcover(:,:,3) - snowcover(:,:,11))./(snowcover(:,:,3) + snowcover(:,:,11));
    % % else
    % %     G = snowcover(:,:,3); G(G==-32768) = NaN;
    % %     SWIR = snowcover(:,:,6); SWIR(SWIR==-32768) = NaN;
    % %     NDSI = (G - SWIR)./(G + SWIR); clear G SWIR;
    % % end
    % % clear snowcover; snowcover = zeros(size(NDSI)); snowcover(NDSI>0.4) = 1; clear NDSI;
    % % 

    %find all IS2 segments for the given date
    inds = find(IS2datenum == IS2datenum_unique(i));
%     IS2track = IS2(IS2datenum == IS2datenum_unique(i)); %index to only data on given date
    IS2X = IS2.Easting(inds); IS2Y = IS2.Northing(inds); 
    
    %find the nearest neighboring snow map pixel & grab its binary snow
    %value (1 = snow, 0 = no snow)
    for j = 1:length(IS2X)
        dist = sqrt((snowXgrid-IS2X(j)).^2 + (snowYgrid-IS2Y(j)).^2); %distance vector
        min_ind = find(dist == min(min(dist)));
        IS2snow(j,1) = snowcover(min_ind);
        clear dist min_ind;
    end

    % %visualize to double-check it worked (can comment-out)
    figure; im = imagesc([snowX./(10^3)],[snowY./(10^3)],snowcover); colormap gray; hold on; %axis xy equal;
    set(gca,'clim',[0 1]);
    plot([IS2X./(10^3)],[IS2Y./(10^3)],'sr','markerfacecolor','none'); hold on;
    plot([IS2X(IS2snow==1)./(10^3)],[IS2Y(IS2snow==1)./(10^3)],'+c'); hold on;
    %geoplot(watershed_border)
    daspect([1 1 1])
    xlabel('Easting [km]')
    ylabel('Northing [km]')
    set(gca,'fontsize',20); set(gca,'Ydir','normal');
    legend('Snow-free ICESat-2 point', 'Snow-covered ICESat-2 point')
    waitfor(im); %close the figure to advance
    
    %add to a vector (to later add to the data table)
    snowvector(inds) = IS2snow;
    snowdatevector(inds) = SnowCoverdates(ix);
    
    %clear out looped variables
    clear d ix snowcover snowX* snowY* inds IS2X IS2Y IS2snow;
end

%% add snow cover flag to table

IS2.snowcover = snowvector;
IS2.snowcoverdate = snowdatevector;
writetable(IS2,[folderpath 'IS2_Data/' outputname]);























