% USAGE: After using bb203, we've got the BP-inferred locations for all
%        qualified aftershocks. Now it's time to calculate slowness errors 
%        and do the calibration for all aftershocks, so as to see how the
%        calibration works for aftershock location.
%        The OUTPUT is same as bb202
% 
% INPUT: 1. mother_loc: Mother directory of all events
%        2. evtlst: event list
%        3. .mat data for runBP (seismic data)
%        4. matched_AFTS_file: matched aftershocks file that has the
%        information of all matched aftershocks: BP-inferred location (La1,Lo1)
%        catalog location (La2,Lo2), Magnitude, Depth, slab interface Depth 
%        at catalog location, origin times in datenum (time2)
%        5. ref: reference station (for calculate slowness errors)


clear
close all
%% PRE-SET

address='/home/liuwei/HPSSD/Supershear/BP2/Workshop2019/';            % address of paths that need to be added below
addpath(genpath(strcat(address,'funcLib/')));   % directory of Back-Projection functions
%%% INPUT
mother_loc = '/home/liuwei/HPSSD/Qinghai/Tele_Cali/TeleDD_Loc/AU/Events';
ap_ca_loc  = '/home/liuwei/HPSSD/Qinghai/Tele_Cali/TeleDD_Loc/AU/Events/ca_ap_loc.mat';
ref = 25;

%     % Parameter set for runBP
%         lon0=119.840; 	% MAINSHOCK lontitude
%         lat0=-0.178;  	% MAINSHOCK latitude
%         dep0 = 10.0;    % MAINSHOCK depth
%         sr=10;          % sample rate
%     	parr=40;      	% start time
%     	begin=0;       	% always 0
%     	over=40;    	% end time
%     	step=1;         % time step
%     	ps=100;          % number of grids for lat
%     	qs=100;          % number of grids for lon
%     	latrange=[-3,1]; % lat range
%     	lonrange=[-1.5,1.5]; % lon range
%     	fl=0.5;        % frequence low of bandpass
%     	fh=2;           % frequence high of bandpass
%     	win=10;         % window length for BP
%         
%         inputband=5;    % i.e. data5.mat
%         Band = 4; % band4: [0.5,2] 
       

cd(mother_loc)
load evtlst.mat
n_evt=length(evtlst_name(:,1));
%% Loop for all events
% for j=1:2
for j = 1 : n_evt
    evtlst_name(j,:)
    cd (evtlst_name(j,:))
    
    load 'Input/Par0.5_2_8.mat'
    load(ap_ca_loc);	% Load ca and ap location of afteshocks
    
    
  	%% STEP 1: Get slowness error for the mainshock.
        
        lat_sta=ret.lat; lon_sta=ret.lon;               % stations
        
        lat_ca = lat_ca; lon_ca = lon_ca; dep_ca(1:9,1) = ret.dep0; % AF catalog
        lat_ap = lat_ap; lon_ap = lon_ap; % AF apparent (BP inferred)
        
        lat_c=mean(lat_ap); lon_c=mean(lon_ap); dep_c =ret.dep0;       % geometrical center
        
        dS2Dplus = get_dS_2Dplus(lat_sta,lon_sta,lat_c,lon_c,dep_c,...
                            lat_ca,lon_ca,dep_ca,lat_ap,lon_ap,ref); % IMPORTANT
        ret.dS0 = dS2Dplus.dS0;
        ret.dSx = dS2Dplus.dSx;
        ret.dSy = dS2Dplus.dSy;
        
 save Input/Par0.5_2_8_cali.mat

    %% runCali
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x0=ret.xori;            % wave data
    ret=rmfield(ret,'xori');% remore ret.xori to save RAM
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
    
    [n,~]=size(x0);         % n: number of stations
    ux=linspace(uxRange(1)+lat0,uxRange(2)+lat0,ps); % lat of grid points
    uy=linspace(uyRange(1)+lon0,uyRange(2)+lon0,qs); % lon of grid points
    
    saveDir=['Input/' dirname '_music_Cali2Dplus_Dir'];
    disp(['Directory: ',saveDir])
    system(['mkdir ' saveDir]);
    cd(saveDir);
    fileID=fopen('logfile','w');
    
    save('parret','ret');
%% Butter Worth Filter   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [BB,AA]=butter(4,[fl fh]/(sr/2));
    for i=1:n
    
        x0(i,:)=filter(BB,AA,x0(i,:));
    end
    for i=1:n 
        x0(i,:)=x0(i,:)/std(x0(i,parr*ret.sr:(parr+30)*ret.sr));
    end
%% Use 1D-ref model to build up time table for grid points
    Ptimesdepth=load('Ptimesdepth.mat');    % load 1D ref-model
    rr=Ptimesdepth.dis;     % distance of 1D model
    dd=Ptimesdepth.dep;     % Depth of 1D model
    tt=Ptimesdepth.Ptt;     % travel times for [rr,dd]

    tlib=zeros(n,ps,qs);    % Initialize the time table
    
   
    for p=1:ps              % Grid point [q,p]
        for q=1:qs          % Grid point [q,p]
            
            lon_cnt=lon_c; % center for this grid
            lat_cnt=lat_c; % center for this grid
            
            dx=111.195*(lat_cnt-ux(p));
            dy=111.195*(lon_cnt-uy(q))*cosd((lat_cnt+ux(p))/2);
            
            sd=phtime2d(lat0,lon0,dep0,ux(p),uy(q),dep0,r(:,2),r(:,1),rr,dd,tt)'...
                        -(ret.dS0 + ret.dSx*dx + ret.dSy*dy); % IMPORTANT;
                    
            tlib(:,p,q)=sd-mean(sd); % tlib is the time table we want
        end
    end
    clear sd
    
%% Generate power field of the source area (Pm) then pick peaks of each sec 
 	
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

        s(i)                              %transfer the serial number back to the frequency
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
            
            disp(['bm bux ' num2str(tmp2(1)) ' buy ' num2str(tmp2(2)) ' max ' num2str(maxp2)]);

         fprintf(fileID,'%f %f %f %f %f %f %f\n',tl,bux,buy,maxp,bux2,buy2,maxp2);

        save ([  num2str(tl-parr-begin) 'smat'],'Pm','Pw');

end

fclose('all');
save('parret','ret');

    % END of current loop
    posteriorNoplotmusic

    close all
    cd ../../..
end