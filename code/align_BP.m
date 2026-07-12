function align_BP(path,Band,ts11,refst,cutoff,plotscale)
% align_BP align seismograms  
%   align_BP(path,Band,ts11,refst,cutoff,plotscale) used after function 
%   'readteleBP' to align seismograms for one earthquake event
%   
%   
%   explaination of the inputs:
%   	path:     	the working address of current project or event, 
%                       e.g. /home/USER/BackProjection/Mexico2017.
%    	Band:           
%     	ts11:     	P arrival time since the straight one (records)
%                       i.e. how long the data start before P-arrival time
%     	refst:   	
%       cutoff:     
%       plotscale:  scaling of the amplitudes of seismograms
%   

%   Copyright 2013-2019 Han Bao & Lingsen Meng's group, UCLA.

filename=strcat(path,'Input/data',num2str(Band-1));
load(filename);
[fl fh win range]=alignband(Band);
ret.x=ret.xori;
load ptimes;
if refst~=0
    lat_ref=ret.r(refst,2);
    lon_ref=ret.r(refst,1);
end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n m]=size(ret.x);

[BB,AA]=butter(4,[fl fh]/(ret.sr/2));
for i=1:n
    ret.x(i,:)=filter(BB,AA,ret.x(i,:));
end

align=seismoAlign();
align.ts11=ts11;% aligned time start from straight one
align.win=win;
align.lt=400;%no change
align.range=range; % upper limitation
align.refst=refst;% first not zero, anyone is ok
align.cutoff=cutoff;
ret=seismoAlign(ret,align);%!!!!!!!important!!!!

ret.scale=plotscale;%no change
ret.lt=300;% no change
ret.ht=350;%no change  above 3 lines for plotting

figure(10);
plotAll1(ret);
figure(11);
plotSta(ret);
if refst~=0
    hold on
    plot( lon_ref, lat_ref, 'r.','MarkerSize',20);
end

filename=strcat(path,'Input/data',num2str(Band));
save(filename,'ret','-v7.3');
end