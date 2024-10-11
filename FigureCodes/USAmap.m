figure(1); clf
usamap("conus");
states = readgeotable('usastatehi.shp');
geoshow(states, 'FaceColor','none')
set(gca,'fontsize',18)