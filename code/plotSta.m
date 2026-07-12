function plotSta(ret)
% function plotSta(ret)
% plot the stations with station name and political maps

if isfield(ret,'stacolor')
    stacolor=ret.stacolor;
else
    stacolor='b.';
end

load haiticoast.dat;
load west_hemi.dat;
load west_political.dat;
load worldcoast.dat
load ptimes;
load japan_coastline.dat;

[nel samlength]=size(ret.x);

r=ret.r;
lat0=ret.opr.lat0;
lon0=ret.opr.lon0;
stanm=ret.nm;


hold on;
% 
% plot(lon,lat,'r*');
plot(lon0,lat0,'r*','MarkerSize',20);
plot( r(:,1), r(:,2),stacolor,'MarkerSize',20);
plot(west_hemi(:,1),west_hemi(:,2),'black');
plot(west_political(:,1),west_political(:,2),'black');
plot(japan_coastline(:,1),japan_coastline(:,2),'black');
plot(worldcoast(:,1),worldcoast(:,2),'black');
axis equal
xlabel('Longitude (^o)');
ylabel('Latitude (^o)');
ylim([-90 90]);
xlim([-180 180]);
box on;
for i=1:nel
end
