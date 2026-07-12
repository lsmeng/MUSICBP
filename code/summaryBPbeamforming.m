% This script will be used to do everything after runteleBP, i.e. movieBP 
% and summaryBP, by that it will pick the peak (radiator) for each second,
% and generate HFdots file. You will need to finish runteleBP or
% runteleBPmusic or runteleBPmusicCali first, then run this script.

% clear all;
% close all;
load parret.mat;

%% load parameter from parret.mat (ret)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lon0=ret.lon0;          % lon of mainshock
    lat0=ret.lat0;          % lat of mainshock
    sr=ret.sr;              % sample rate
    parr=ret.parr;          % P-arrival time
    tend=ret.end-ret.begin; % duration of movie
    step=ret.step;          % time step
    ps=ret.ps;              % number of lat grids
    qs=ret.qs;              % number of lon grids
    uxRange=ret.latrange;   % range of latitude
    uyRange=ret.lonrange;   % range of longitude
%% Some other pre-set  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ncontour=200;           
    load worldcoast.dat;
    load asiapolitical;
    t=1:step:tend;          % time
    Pm=zeros(ps,qs);        % MUSIC power
    ux=linspace(uxRange(1)+lat0,uxRange(2)+lat0,ps); % lat of grid points in source area   
    uy=linspace(uyRange(1)+lon0,uyRange(2)+lon0,qs); % lon of grid points in source area
    bux=zeros(length(t),10);                         % bux will be the lat of peaks
    buy=zeros(length(t),10);                         % buy will be the lon of peaks
    Power=zeros(length(t),10);                       % Power will be the beamforing power (Pm) of peaks
    nw=zeros(1,length(t));                           % nw will be the number of peaks

    x1=zeros(length(t),2);          % x1 will be the lat of the man_pick peak(s)
    y1=zeros(length(t),2);          % y1 will be the lon of the man_pick peak(s)

%% Pre-process of 0smat.mat (build up the loop of GIF, and plot 0s movie)
    h=figure(2);
    set(gcf,'Position',[100 1000 200*(uyRange(2)-uyRange(1)) 200*(uxRange(2)-uxRange(1))]);
    hold on;
    
    %%%%%%%%%%%%%% plot 0s movie and build up the loop of GIF %%%%%%%%%%%%%
    load('0smat.mat');
    [c,ch]=contourf(uy,ux,real(Pm)/20,ncontour); % Plot the MUSIC power (Pm)
    axis equal;
    set(ch,'edgecolor','none');
    colormap(jet);
    colorbar;
    xlabel('^oE');
	ylabel('^oN');
    chsize(20) % set FontSize of Label and Title
    drawnow;

    f = getframe(h);
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    clf;    % clean current figure. 
    imwrite(im,map,'movie.gif','DelayTime',0.50,'LoopCount',inf)
    kt=0;   % Name-value pair 'LoopCount',Inf causes the animation to continuously loop.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Process every second after 0s
for j=1:length(t)
    kt=kt+1;
    
    hold on;
    title([num2str(t(j)) 's']);
    load( [num2str(t(j)) 'smat.mat']);  % loat 't'smat.mat
    
    Pm=abs(real(Pm));
    Pm=Pm-min(Pm(:));   % MUSIC Power (Pm)

	%% Pick 'localMaximum' than use 'sortrows' to sort peaks with decreasing Pm
        [pw,qw]=localMaximum(abs(Pm),[2 2]); % pw,qw: index of lon/lat of peaks
        clear tmp;
        for jj=1:length(pw)
            tmp(jj)=abs(Pm(pw(jj),qw(jj)));
        end
        %%%%%%%%%%% Sort peaks with the deceasing energy (tmp) %%%%%%%%%%%%
        ttmp=sortrows([pw qw tmp'],-3); % 
        pw=ttmp(:,1);
        qw=ttmp(:,2);
        tmp=ttmp(:,3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%% DETERMINE if peak(jj) is too small (<0.5*tmp(1)) %%%%%%%%%
        nw(kt)=min([ 3 length(tmp) ]);% nw will be the number of peaks 
        for jj=1:nw(kt)         
            if tmp(jj)<0.5*tmp(1)   % note that tmp(1) is the max power at this time
                nw(kt)=jj-1;        % if so, than delete it, as we decrease nw.
                break;              % nw is the number of peaks
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%% Value Assignment of lat/lon of peaks %%%%%%%%%%%%%%%
        bux(kt,1:nw(kt)) = interp1(1:length(ux),ux,pw(1:nw(kt)));
        buy(kt,1:nw(kt)) = interp1(1:length(uy),uy,qw(1:nw(kt)));
        Power(kt,1:nw(kt))=interp2(uy,ux,abs(Pm),buy(kt,1:nw(kt)),bux(kt,1:nw(kt)));% power of the peaks
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot Energy distribution and autopick_peaks (BP image)

        [c,ch]=contourf(uy,ux,real(Pm)/20,ncontour);    % Plot the MUSIC power (Pm)
        axis equal;
        set(ch,'edgecolor','none');
        colormap(jet);
        colorbar;
        xlabel('^oE');
        ylabel('^oN');
        chsize(20) % set FontSize of Label and Title
        drawnow;
        
        plot(lon0,lat0,'r*','MarkerSize',10);   % plot the Epicenter
        text(lon0,lat0,'\leftarrow Epicenter','Color','red');
        for jj=1:nw(kt)
%             text(buy(kt,jj),bux(kt,jj),num2str(jj),'Color','cyan','HorizontalAlignment','center')
        end
        
    %% Save every second images into a GIF file. %%%%%%%%%%%%%%%%%%%%%%%%%%
    f=getframe(h);
    im= rgb2ind(f.cdata,map,'nodither');
    clf
    imwrite(im,map,'movie.gif','DelayTime',0.5,'WriteMode','append') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
save movieBP;

%% Write HFdots file
Power = sqrt(Power); % the square root of amplitude
Power = Power / max(Power(:));

load ptimes;
for j=1:length(t)
    sd=phtime(1,lat0,lon0,bux(j,1),buy(j,1),ret.r(:,2),ret.r(:,1),rr,tt)'; 
    t(j)=t(j)-mean(sd);
end

x3=[ t' bux(1:length(t),1) buy(1:length(t),1) Power(1:length(t),1)];
save('HFdots','x3','-ascii'); % generate the HFdots file

%% SummaryteleBP

dep0=ret.dep0;
h12=figure(2);
set(h12, 'Position', [100, 100, 600, 600]);
set(gcf, 'PaperPositionMode', 'auto') ;
cmm=colormap;
[ nmm, tmp0]=size(cmm);
hold on;
c0=linspace(1,nmm,tend);
r=interp1(1:nmm,cmm(:,1),c0);
g=interp1(1:nmm,cmm(:,2),c0);
b=interp1(1:nmm,cmm(:,3),c0);
r1=interp1(1:tend,r,t);
g1=interp1(1:tend,g,t);
b1=interp1(1:tend,b,t);

kt=0;
for j=1:length(t)
    kt=kt+1;
    for jj=1:1
        marker='o';
        scatter(buy(kt,jj),bux(kt,jj),abs(Power(kt,jj))*100,[r1(j) g1(j) b1(j)],marker,'filled');
    end
end
plot(worldcoast(:,1),worldcoast(:,2),'black','LineWidth',2);
%%% Plot Palu Fault
    load Palu_fault.dat
    plot(Palu_fault(:,1),Palu_fault(:,2),'red','LineWidth',2);
hold on
    JPtrench=[
  142.0400   34.5000
  142.1000   35.0000
  142.2000   35.3400
  142.4200   35.6700
  142.8500   36.1900
  143.2200   36.6500
  143.5100   37.1400
  143.7900   37.6600
  143.9500   38.1600
  144.1600   39.0600
  144.2900   39.6900
  144.4100   40.2500
  144.5000   40.6700
  144.7000   41.1000
   ] 
   plot(JPtrench(:,1),JPtrench(:,2),'b-.','LineWidth',2);
 
    
    
    plot(ret.lon0,ret.lat0,'r*','MarkerSize',10);
set(gca,'DataAspectRatio',[1/cosd(lat0) 1 1])
ylim([min(ux) max(ux)]);
xlim([min(uy) max(uy)]);
colorbar;
colormap(cmm);
caxis([0 tend]);
xlabel('Longitude (^o)');
ylabel('Latitude (^o)');
title('Summary');
box on;
print('-dpdf','-r300','summary_BP.pdf'); 

    %% Plot Power-time
    figure;plot(real(Power(:,1)) );
    set(gcf, 'Position', [100, 100, 600, 600]);
    set(gcf, 'PaperPositionMode', 'auto') ;
    box on; 
    xlabel('Time (s)');
    ylabel('Normalized Power');
    title('Power vs time');
    print('-dpdf','-r300','power_BP.pdf'); 

    %% Plot distance vs time projection along strike
    figure;
epi_dis=distance22(bux(3,1),buy(3,1),bux(:,1),buy(:,1))*110.49;
epi_az=azimuth(bux(3,1),buy(3,1),bux(:,1),buy(:,1));
pjt_az=0;
pjt_dis=epi_dis.*cosd(epi_az-pjt_az);
pjt_pen_dis=-epi_dis.*sind(epi_az-pjt_az);
% x3=[ t' bux(1:length(t),1) buy(1:length(t),1) Power(1:length(t),1)/max(Power(1:length(t),1)) epi_dis(1:length(t)) epi_az(1:length(t)) pjt_dis(1:length(t)) pjt_pen_dis(1:length(t))];
% 
% save('HFdots','x3','-ascii');

%h12=figure(13);
hold on;
set(gcf, 'Position', [100, 100, 600, 600]);
set(gcf, 'PaperPositionMode', 'auto') ;

% az0=40;az1=160;
% N=180;
% c0=linspace(1,64,N);
% r=interp1(1:64,cmm(:,1),c0);
% g=interp1(1:64,cmm(:,2),c0);
% b=interp1(1:64,cmm(:,3),c0);
% r1=interp1(linspace(az0,az1,N),r,epi_az);
% g1=interp1(linspace(az0,az1,N),g,epi_az);
% b1=interp1(linspace(az0,az1,N),b,epi_az);
for j=1:length(t)
    kt=kt+1;
    for jj=1:1
        marker='o';
        scatter(-pjt_dis(j),t(j),abs(Power(j,jj))*150,[r1(j) g1(j) b1(j)],marker,'filled');
    end
end
colorbar;
colormap(cmm);
caxis([0 tend]);
xlabel('Distance (km)');
ylabel('Time (s)');
title('Distance along N-S direction')
box on;
print('-dpdf','-r300','dis_NS_BP.pdf'); 

figure;
hold on;
set(gcf, 'Position', [100, 100, 600, 600]);
set(gcf, 'PaperPositionMode', 'auto') ;
for j=1:length(t)
    kt=kt+1;
    for jj=1:1
        marker='o';
        scatter(epi_dis(j),t(j),abs(Power(j,jj))*150,[r1(j) g1(j) b1(j)],marker,'filled');
    end
end
colorbar;
colormap(cmm);
caxis([0 tend]);
xlabel('Distance (km)');
ylabel('Time (s)');
title('Distance to mainshock epicenter')
box on;
print('-dpdf','-r300','dis_epi_BP.pdf');
