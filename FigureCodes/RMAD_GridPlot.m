%% Inputs
clearvars;
addpath(['/Users/karinazikan/Documents/ICESat2-AlpineSnow/functions'])
addpath(['/Users/karinazikan/Documents/cmocean'])

%site abbreviation for file names
abbrev = 'RCEW';
%Product abbreviation for files
prod_abbrev = 'A6-40';
%Folder path
folderpath = ['/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/' abbrev '/'];

%Turn slope correction off or on
slope_correction = 0; % 0 = off, 1 = on

%set colors
colors{1} = cmocean('-dense',6);
colors{2} = cmocean('-algae',5);
colors{3} = cmocean('ice',5);
colors{4} = cmocean('-amp',5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
% IS2 data
filepath = [folderpath 'IS2_Data/' prod_abbrev '/RCEW_A6-40-20210123_rmadGrid.csv'];
rmad_grid = readmatrix(filepath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
A1 = [-8:8];
figure(1); clf; hold on
im = imagesc(rmad_grid);
xline(9); yline(9);
ax = gca;
ax.DataAspectRatio = [1 1 1];
set(gca,'fontsize',18,'xlim',[0.5,17.5],'ylim',[0.5,17.5])
xticks(1:length(A1)); yticks(1:length(A1));
xticklabels(A1); yticklabels(A1);
xlabel('Easting offset (m)'); ylabel('Northing Offset (m)'); 
c = colorbar;
c.Label.String = 'NMAD (m)';