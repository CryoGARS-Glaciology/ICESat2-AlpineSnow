clear all

%Folder path 
folderpath = '/Users/karinazikan/Documents/GitHub/ICESat2-AlpineSnow/Sites/RCEW/';
DTM_name = 'RCEW_1m_WGS84UTM11_WGS84.tif';
%site abbreviation for file names
abbrev = 'RCEW';
%Read in files
file = readtable([folderpath 'IS2_Data/' abbrev '-ICESat2-ATL03']);
[DTM,Ref] = readgeoraster([folderpath 'DEMs/' DTM_name]);

%DEMgrid
if isfield(Ref,'LatitudeLimits')
    [latgrid,longrid] = meshgrid(Ref.LongitudeLimits(1)+0.5*Ref.CellExtentInLongitude:Ref.CellExtentInLongitude:Ref.LongitudeLimits(2)-0.5*Ref.CellExtentInLongitude,...
        Ref.LatitudeLimits(2)-0.5*Ref.CellExtentInLatitude:-Ref.CellExtentInLatitude:Ref.LatitudeLimits(1)+0.5*Ref.CellExtentInLatitude);
    [xgrid, ygrid,~] = wgs2utm(latgrid,longrid);
else
    x = Ref.XWorldLimits(1)+0.5*Ref.CellExtentInWorldX:Ref.CellExtentInWorldX:Ref.XWorldLimits(2)-0.5*Ref.CellExtentInWorldX;
    if strcmp(Ref.ColumnsStartFrom,'north')
        y = Ref.YWorldLimits(2)-0.5*Ref.CellExtentInWorldY:-Ref.CellExtentInWorldY:Ref.YWorldLimits(1)+0.5*Ref.CellExtentInWorldY;
    else
        y = Ref.YWorldLimits(1)+0.5*Ref.CellExtentInWorldY:Ref.CellExtentInWorldY:Ref.YWorldLimits(2)-0.5*Ref.CellExtentInWorldY;
    end
    [xgrid, ygrid] = meshgrid(x, y); % create grids of each of the x and y coords
end
DTM(DTM>(2.5*10^3)) = NaN;


%% visualize data
figure(1);
imagesc(x,y,DTM);
hold on
scatter3(file,"X","Y","Z",Marker=".");

% %% Identify tracks that dont intersect
% dates = datetime(file.time.Year,file.time.Month,file.time.Day);
% dates_unique = unique(dates);
% for i = 1:length(dates_unique)
%     ix = find(dates == dates_unique(i));
%     tracks{i} = file(ix,:);
% end
% % for i = 1:(length(tracks)-1)
% %     for k = (i+1):length(tracks)
% %         check{i,k} = polyxpoly(tracks{i}.X,tracks{i}.Y,tracks{k}.X,tracks{k}.Y,'unique');
% %     end
% % end
% 
% figure(2);
% hold on
% for i = 1:length(tracks)
% scatter3(tracks{i},"X","Y","time",Marker=".");
% end
% legend(datestr(dates_unique),'FontSize',12)
% 
% %% Remove tracks that dont intersect
% ix = find(dates ~= '2019-08-26');
% file2 = file(ix,:);
% 
% figure(3);
% scatter3(file2,"X","Y","Z",Marker=".");
% 
% %% Write new point cloud
% writetable(file2,[folderpath 'IS2_Data/' abbrev '-ICESat2-ATL03-2.csv']);

%% Make 11m elevation array
x11 = downsample(x,11); y11 = downsample(y,11);

for i = 1:length(x11)
    for j = 1:length(y11)
        ix = find(file.X <= (x11(i)+ 5.5 & file.X >= (x11(i)-5.5) & file.Y <= (y11(i)+ 5.5 & file.Y >= (y11(i)-5.5)))); 
        if isempty(ix) == 1
            IS2grid(j,i) = NaN;
        else
            IS2grid(j,i) = min(file.Z(ix));
        end
    end
end

%% Create the referencing matrix (R)
R = maprasterref('RasterSize',size(IS2grid), ...
    'XLimWorld',[min(x11) max(x11)],'YLimWorld',[min(y11) max(y11)]);
R.ProjectedCRS  = Ref.ProjectedCRS; 
R.ColumnsStartFrom = 'north';

%% Export raster
filename = [folderpath 'DEMs/' abbrev '-ICESat2-ATL03.tif'];
geotiffwrite(filename,IS2grid,R,'CoordRefSysCode','epsg:32611')