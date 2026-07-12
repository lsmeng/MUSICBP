
function readteleBP(path,lon0,lat0,sr,ori,display,Preshift,plotscale)
%readteleBP read .SAC files of an earthquake event.
%   readteleBP(path,lon0,lat0,sr,display,Preshift) read .SAC files that
%   restored in path/Data. 'path' is the working address of current project
%   or event, e.g. /home/USER/BackProjection/Mexico2017.
%   
%   explaination of the inputs:
%   	path:           the working address of current project or event, 
%                         e.g. /home/USER/BackProjection/Mexico2017.
%    	(lon0,lat0):    the epicenter location of the current event
%     	sr:             is sample rate.
%       
%     	display:        time length of waveform showing by plotAll1(ret)
%       plotscale:    	scaling of the amplitudes of seismograms
% 
%
%   Copyright 2013-2017 Han Bao & Lingsen Meng's group, UCLA.


opr=readAllSac();
opr.bpbool=false;
opr.lon0=lon0;
opr.lat0=lat0;
opr.bp=[0.01 2];
opr.snrFilterbool=false;
opr.snrFilter=[0.1 0.1 2 -20 -10 100 130];
opr.sr=sr; 
opr.ori=ori;
datalog=strcat(path,'Data/');
ret=readAllSac(datalog,opr);
%
ret.xori=ret.x;
%save('data','ret','-v7.3');

load ptimes;
[n m]=size(ret.x);
fprintf ('The number of stations is %d\n',n);
%I=2:8;
figure(1);
plotSta(ret);
%
load ptimes
if Preshift==true
    shift=interp1(rr,tt,ret.rdis);
    for i=1:n
        ret.x(i,:)=specshift(ret.x(i,:),ret.sr*(shift(i)-shift(1)));
    end
end
ret.xori=ret.x;
fl=0.1;
fh=0.5;
[BB,AA]=butter(4,[fl fh]/(ret.sr/2));
for i=1:n
    ret.x(i,:)=filter(BB,AA,ret.x(i,:));
end
ret.x=ret.x(:,1:display*ret.sr);
for i=1:n
    ret.x(i,:)=ret.x(i,:)/std(ret.x(i,10*ret.sr:display*ret.sr));
    ret.xori(i,:)=ret.xori(i,:)/std(ret.xori(i,10*ret.sr:display*ret.sr));
end
ret.scale=plotscale;
ret.lt=300;
ret.ht=350;
% [delete_num delete_i]=size(ret.nm);
% delete_c=[19:55,1];
% B=delete_a(delete_num,delete_c');
% ret=orderAll(ret,0,B);
I = find(ret.rdis<90);
ret=orderAll(ret,0,I);
figure(10);
plotAll1(ret);
filename=strcat(path,'Input/data0');
save(filename,'ret','-v7.3');
end