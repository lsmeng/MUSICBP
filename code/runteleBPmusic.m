function runteleBPmusic(path,BPfile)
% If had bug, add path to MUSICBP
addpath(pwd)
load(BPfile)
load ptimes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=ret.xori;            % wave data
r=ret.r;                % stations' location
lon0=ret.lon0;          % hypocenter longitude
lat0=ret.lat0;          % hypocenter latitude
dep0=ret.dep0;          % hypocenter depth
sr=ret.sr;              % sample rate
parr=ret.parr;          % start time
begin=ret.begin;        % always 0
over=ret.end;           % end time
step=ret.step;          % time step
ps=ret.ps;              % number of grids for lat
qs=ret.qs;              % number of grids for lon
uxRange=ret.latrange;   % lat range
uyRange=ret.lonrange;   % lon range
fl=ret.fl;              % frequence low of bandpass
fh=ret.fh;              % frequence high of bandpass
win=ret.win;            % window length for BP
dirname=ret.dirname;    % name of directory that will be created
Nw=ret.Nw;
fs=ret.fs;
%%%%%%%%%%%%%%%%%%%%%%%

[n m]=size(x0);
ux=linspace(uxRange(1)+lat0,uxRange(2)+lat0,ps); 
uy=linspace(uyRange(1)+lon0,uyRange(2)+lon0,qs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd([path 'Input']);
fprintf dirname
saveDir=[dirname '_MUSIC_Dir'];
system(['mkdir ' saveDir]);
cd(saveDir);
fileID=fopen('logfile','w');

rmfield(ret,'xori');
save('parret','ret');
%% Butter Worth Filter  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
[BB,AA]=butter(4,[fl fh]/(sr/2));
for i=1:n
    
    x0(i,:)=filter(BB,AA,x0(i,:));
end
for i=1:n 
        x0(i,:)=x0(i,:)/std(x0(i,parr*ret.sr:(parr+30)*ret.sr));
end
%% Use 1D-ref model to build up time table for grid points  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(ret,'pP')
	Ptimesdepth=load('pPtimesdepth.mat'); % load 1D ref-model
else
Ptimesdepth=load('Ptimesdepth.mat');      % load 1D ref-model
end
rr=Ptimesdepth.dis;                       % distance of 1D model
dd=Ptimesdepth.dep;                       % Depth of 1D model
tt=Ptimesdepth.Ptt;                       % travel times for [rr,dd]
tlib=zeros(n,ps,qs);                      % Initialize the time table

for p=1:ps                                % Grid point [q,p]
    
    for q=1:qs                            % Grid point [q,p]
        
              sd=phtime2d(lat0,lon0,dep0,ux(p),uy(q),dep0,r(:,2),r(:,1),rr,dd,tt)';% %Calculate time difference
   
        tlib(:,p,q)=sd-mean(sd);          % tlib is the time table we want
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%

for tl=parr+begin:step:parr+over          %for each second
    
    
    th=tl+win;                            %in time window
    display(['t=' num2str(tl-parr-begin) 's']); %show the time


    Pm=zeros(ps,qs);                      %initialize the power space
    Pw=zeros(ps,qs);                      %initialize the power space

    S=cmtmall(x0(:,tl*sr:th*sr-1),Nw);    %S is a n*n autocorrelation matrix, contain the waveform information received by stations
    s=linspace(0,sr-sr/((th-tl)*sr),(th-tl)*sr);%length(s)=win*sr, s is used for the loop
    fli=round(interp1(s,1:length(s),fl)); %the serial number for low limit frequency
    fhi=round(interp1(s,1:length(s),fh)); %the serial number for high limit frequency
%%  for every frequency in the range of fl to fh
    for i=fli:fs:fhi

        s(i);                              %transfer the serial number back to the frequency
        Pm1=zeros(ps,qs);                 %Pm1 is the part of Pm in this frequency
        Pw1=zeros(ps,qs);                 %Pw1 is the part of Pw in this frequency
        clear Uv A Un a wi;
        Rxx=zeros(n,n);
        for j=1:n
            for k=1:n
                Rxx(j,k)=S(i,j,k);        %get the autocorrelation matrix
            end
        end
        [Uv,A]=eig(Rxx);                  %get the eigenvectors
        As=zeros(n,n);
        un=zeros(n,n);
        us=zeros(n,n);
        M=rank(Rxx);
        
        
        
        un(:,1:n-M)=Uv(:,1:n-M);          %un is the useful signal
        
        Un=un*un';
        for p=1:ps
            
            for q=1:qs
                
                a=exp(-1i*2*pi*s(i)*tlib(:,p,q));% a is a n*1 vector
                Pm1(p,q)=((a'*a)/(a'*Un*a));     
                Pw1(p,q)=((a'*Rxx*a)/(a'*a));
            end
        end
        
        Pm=Pm+Pm1;    %sum the Pm1 of in the frequency band 
        Pw=Pw+Pw1;    %sum the Pw1 of in the frequency band 
        
    end
 %%%%%%%%  Pick power peak for the power field at tl(s) %%%%%%%%%%%%
           tmp1=peakfit2d(real(Pm));
            bux=interp1(1:length(ux),ux,tmp1(1));       % lat of peak at time tl
            buy=interp1(1:length(uy),uy,tmp1(2));       % lon of peak at time tl
            maxp=interp2(1:ps,1:qs,Pm',tmp1(1),tmp1(2),'linear',0);% power of peak
            
             tmp2=peakfit2d(real(Pw));
            bux2=interp1(1:length(ux),ux,tmp2(1));      % lat of peak at time tl
            buy2=interp1(1:length(uy),uy,tmp2(2));      % lon of peak at time tl
            maxp2=interp2(1:ps,1:qs,Pw',tmp2(1),tmp2(2),'linear',0);% power of peak
            
           % disp(['bm bux ' num2str(tmp2(1)) ' buy ' num2str(tmp2(2)) ' max ' num2str(maxp2)]);

         fprintf(fileID,'%f %f %f %f %f %f %f\n',tl,bux,buy,maxp,bux2,buy2,maxp2);

        save ([  num2str(tl-parr-begin) 'smat'],'Pm','Pw');

end
rmfield(ret,'xori');
save('parret','ret');
end