function runteleBPbmfm(path,BPfile)
% If had bug, add path to MUSICBP
addpath(pwd)
load(BPfile)
load ptimes;
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

    [n,~]=size(x0);         % n: number of stations
    ux=linspace(uxRange(1)+lat0,uxRange(2)+lat0,ps); % lat of grid points
    uy=linspace(uyRange(1)+lon0,uyRange(2)+lon0,qs); % lon of grid points
    
    cd([path 'Input']);
    fprintf dirname
    saveDir=[dirname '_bmfm_Dir'];
    system(['mkdir ' saveDir]);
    cd(saveDir);
    fileID=fopen('logfile','w');
    
    save('parret','ret');
%% Butter Worth Filter   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [BB,AA]=butter(4,[fl fh]/(sr/2));
    for i=1:n
        x0(i,:)=filter(BB,AA,x0(i,:));
    end
    
%% Use 1D-ref model to build up time table for grid points
    Ptimesdepth=load('Ptimesdepth.mat');    % load 1D ref-model
    rr=Ptimesdepth.dis;     % distance of 1D model
    dd=Ptimesdepth.dep;     % Depth of 1D model
    tt=Ptimesdepth.Ptt;     % travel times for [rr,dd]

    tlib=zeros(n,ps,qs);    % Initialize the time table
    
    for p=1:ps              % Grid point [q,p]
        for q=1:qs          % Grid point [q,p]
            sd=phtime2d(lat0,lon0,dep0,ux(p),uy(q),dep0,r(:,2),r(:,1),rr,dd,tt)';
            
            tlib(:,p,q)=sd-mean(sd); % tlib is the time table we want
            if max(abs(sd)) > 200
                disp('ERROR: rdis is too long')
                continue
            end
        end
    end
    clear sd
    
%% Generate power field of the source area (Pm) then pick peaks of each sec 
 	Pm=zeros(ps,qs);        % Initialize the power table for grid points
  	
    for tl=parr+begin:step:parr+over    % for each time window [tl,th]
        th=tl+win;                      % th is the end of each window
        disp(['t=' num2str(tl-parr-begin) 's']); % display current second
        y=zeros(n,win*sr);
        
        %%%%%%%%  Using Sliding Window Beamforming to get energy  %%%%%%%%%
        for p=1:ps       	% Grid point [q,p]
            for q=1:qs   	% Grid point [q,p]

                sd=tlib(:,p,q)';  % time table of [q,p] to all stations
                for k=1:n	% !!! y(k,:) are all the samples from tl to th
                            % for trace(station k), assuming [q,p] is the
                            % source point.
                    y(k,:)=x0( k, floor((tl+sd(k))*sr) : floor((th+sd(k))*sr-1) );
                end

                Pm(p,q)=sum(sum(y(:,:),1).^2); 
                            % one value for point [q,p] at time tl: summation
                            % of n*(win*sr), i.e. all samples in the time
                            % window of all the n stations
            end
        end
        %%%%%%%%  Pick power peak for the power field at tl(s) %%%%%%%%%%%%
        tmp1=peakfit2d(real(Pm));
        bux=interp1(1:length(ux),ux,tmp1(1)); % lat of peak at time tl
        buy=interp1(1:length(uy),uy,tmp1(2)); % lon of peak at time tl
        maxp=interp2(1:ps,1:qs,Pm',tmp1(1),tmp1(2),'linear',0); % power of peak
           
        %disp(['bm bux ' num2str(tmp1(1)) ' buy ' num2str(tmp1(2)) ' max ' num2str(maxp)]);

        fprintf(fileID,'%f %f %f %f\n',tl,bux,buy,maxp); % write into 'logfile'
        save ([  num2str(tl-parr-begin) 'smat'],'Pm');
    end

fclose(fileID);
save('parret','ret','-v7.3')