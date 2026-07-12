function t=phtime(f,lat1,lon1,lat2,lon2,lats,lons,rr,tt)
% f=1;
% lat1=lat0;
% lon1=lon0;
% lat2=ux(p);
% lon2=uy(q);
% tt=ttime;
% lats=lat;
% lons=lon;
r1=distance22(lat1,lon1,lats,lons);

t1=interp1(rr,tt,r1);

r2=distance22(lat2,lon2,lats,lons);

t2=interp1(rr,tt,r2);

t=t2-t1;

%ph=exp(-1i*2*pi*f*t);
%if max(r1)>90||max(r2)
% r=0:0.1:30;
% t=zeros(1,length(r));
% for i=1:length(r)
% rr=num2str(r(i));
% [status,result]=system(['phtimes 10 ' rr ' P']);
% t(i)=str2double(result(5:20));
% endfigu