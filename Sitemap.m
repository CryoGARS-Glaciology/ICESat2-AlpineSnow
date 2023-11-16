%%% Make sitemap figure
%%%
%%% SPECIFIED INPUTS:
%%%     site = path to watershed outline
%%%     
%%% OUTPUTS:
%%%     
%%%
%%% Last updated: Nov 2023 by Karina Zikan

%% Inputs
site = readgeotable("/Users/karinazikan/Documents/ICESat2-AlpineSnow/Sites/RCEW/ROIs/RCEW-outline_WGS84.shp");

%%
states = readgeotable('usastatehi.shp');
row = states.Name == "Idaho";
state = states(row,:);

%%
figure(1); clf
geoplot(site,'linewidth',2,'Color','r')
geobasemap colorterrain
hold on
geoplot(state,'FaceColor','none','linewidth',2)
set(gca,'fontsize',16);
%legend('Reynold Creek Experimental Watershed')
