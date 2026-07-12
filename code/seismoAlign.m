function ret = seismoAlign(ret,align)
% function ret = seismoAlign(ret,align)
% function align the data matrix ret.x with first arrivals by
% crosscorrelate all the seismogram with a reference seismogram

def=struct('win',9,...% length of the window for the cross correlation
    'lt',100,...% range*2 is devided by lt as the step of the cross correlation
    'ts11',140,...% starting time of the window
    'range',10,...% plus and minus range of the possible time shift
    'cutoff',0,...% a correlation coefficient cut off below which will be throw away
    'refst',1,...% the number of the reference seismograms
    'bpbool',0,...%whether to do the cross correlation within a frequency band
    'bp',[0.5 2],...% the frequency band
     'timeShiftKnown',0,...,%input know timeshift default 0:unknown timeshift. a vector known timeshift suppress the timeshift caculated.
     'shiftT1bool',false);% use theoretical arrival time to align data
                         % can be get by gawk '{print "saclst t1 f 1_*"$1}' filelist1 | sh > arrival;
if ~nargin %user wants default options, give it and stop
    ret= def;
    return
elseif nargin<2
    align=def;
end


[nel lengthSam]=size(ret.x);
x=ret.x;
x0=ret.x;
sr=ret.sr;
[BB,AA]=butter(4,align.bp/(sr/2));

if align.bpbool==true
    for i=1:nel
        x(i,:)=filter(BB,AA,x(i,:));
    end
   
end

%parameter settings

win=align.win;%9
lt=align.lt;%200
ts11=align.ts11;%ori
range=align.range;%10
cutoff=align.cutoff;%0.85
refst=align.refst;%1
  timeshift=zeros(1,nel);
   t3=(0:length(x)-1)/sr;
 timeShiftKnown=align.timeShiftKnown;
 shiftT1bool=align.shiftT1bool;
 B=1:nel;
 
 if timeShiftKnown(1)~=0% unknowtimeshift need to caculate
     
     timeshift=timeShiftKnown;
 elseif shiftT1bool==true;
     timeshift=ret.t1-200;
 else
     if refst==0
         xx0=mean(x(:,ts11*sr:(ts11+win)*sr),1);
     else
         xx0=x(refst,ts11*sr:(ts11+win)*sr);%reference seismograms
     end
     tou=linspace(-range,range,lt);% possible time shift range
     
     xcr=zeros(nel,lt);% all correlation coeff
     mxcr=zeros(1,nel);% maxium cross corelation coeff
     
     for i=1:nel
         i
         for j=1:lt
             clear tip xx;
             ts22=ts11+tou(j);
             tip=linspace(ts22,ts22+win,win*sr+1);
             xx=interp1(t3,x(i,:),tip);
             xcr(i,j)=sum(xx0.*xx)/sqrt((sum(xx0.^2)*sum(xx.^2)));%-sum((xx0-xx).^2)/sum((xx0).^2);
         end
     end
     
     
     for i=1:nel
         mma=max(xcr(i,:));
         if mma<cutoff&&~isnan(cutoff)
             B(i)=0;
             disp('cut');
         end
         for j=1:lt
             if abs(xcr(i,j))==mma
                 timeshift(i)=tou(j);
                 mxcr(i)=xcr(i,j);
             end
         end
     end
     
 end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     shifting the data with y2;
timeshift=timeshift-timeshift(1);
    x1=zeros(nel,length(x0));
    for i=1:nel
           x1(i,:)=specshift(x0(i,:),timeshift(i)*sr);
%                      x1(i,:)=mxcr(i)/abs(mxcr(i))*specshift(x0(i,:),timeshift(i)*sr);
    end
    if isfield(ret,'xori')
        for i=1:nel
            ret.xori(i,:)=specshift(ret.xori(i,:),timeshift(i)*sr);
            %                      x1(i,:)=mxcr(i)/abs(mxcr(i))*specshift(x0(i,:),timeshift(i)*sr);

            %disp(['timeshift' num2str(timeshift(i))]);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%
    x0=x1;%change both
    %%%%%%%%%%%%
        ret.timeshift=ret.timeshift+timeshift;% the timeshift 

    B=find(B(:)~=0);

      ret=orderAll(ret,0,B);

     ret.x=x0(B,:);
%       if isfield(ret,'xori')
%      ret.xori=ret.xori(B,:);
%       end
%     ret.r=ret.r(B,:);
%     ret.nm=ret.nm(B,:);
%        if isfield(ret,'rdis')
%      ret.rdis=ret.rdis(B);
%       end
% %     ret.az=ret.az(B);
%        if isfield(ret,'x1')
%     ret.x1=x1(B,:);
%        end
%     ret.t1=ret.t1(B);% the aligned the seismograms
    if exist('mxcr','var')
    ret.mxcr1=mxcr;
    end
    %%%%%%%%%%%%%%%%%%%%%%
    
 
    
    
    
    