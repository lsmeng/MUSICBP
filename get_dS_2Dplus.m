% USAGE: this function will get the 2D slowness error with G constant term,
%        i.e. dS0, dSx, and dSy.
% INPUT: 1. lat_sta,lon_sta: stations' location
%        2. lat0,lon0,dep0: mainshock location and depth
%        3. lat_ca,lon_ca,dep_ca: aftershocks location and depth (catalog)
%        4. lat_ap,lon_ap: apparent aftershocks location and depth (apparent BP inferred)
%
% OUTPUT: ds2Dplus (G structure): having stations location and
%         corresponding slowness errors (dS0, dSx, and dSy)
%
% Brandon Han Bao (hbrandon@ucla.edu)  Aug 29th, 2018
% Liuwei Xu, add line 54 "dS0=dS0-mean(dS0);", Nov 25th, 2021

function dS2Dplus=get_dS_2Dplus(lat_sta,lon_sta,lat0,lon0,dep0,lat_ca,lon_ca,dep_ca,lat_ap,lon_ap,ref)

n_sta = length(lat_sta); 	% number of stations
n_evt = length(lat_ca);     % number of aftershocks

%% Calculate dt (differential time) between apparent BP location and catalog location of aftershock(i) to all stations
    Ptimesdepth=load('Ptimesdepth.mat');
    rr=Ptimesdepth.dis;
    dd=Ptimesdepth.dep;
    tt=Ptimesdepth.Ptt;

    dx=zeros(n_evt);dy=dx;sd1=zeros(n_evt,n_sta);sd2=sd1;sd=sd1; % initialization
    for i=1:n_evt
        %% For 2D, no need azimuth
        dx(i)= 111.195*(lat0-lat_ap(i));  % dx(i) is the catalog North-Sorth-distance btwn main and aftshk(i)
      	dy(i)= 111.195*(lon0-lon_ap(i))*cosd((lat0+lat_ap(i))/2);
        sd1(i,:)=phtime2d(lat0,lon0,dep0,lat_ap(i),lon_ap(i),dep_ca(i),lat_sta,lon_sta,rr,dd,tt)';% travel time from P along ray path
        sd2(i,:)=phtime2d(lat0,lon0,dep0,lat_ca(i),lon_ca(i),dep_ca(i),lat_sta,lon_sta,rr,dd,tt)';% P along ray path
        sd(i,:)=sd2(i,:)-sd1(i,:);
    end

%% Build the Matrix D[ddt(n_evt*n_sta,1)] and G[dX(n_evt*n_sta,n_sta*3)].   
%%% where M[ds(n_sta*3,1)] = D\G, or D=G*M.     M is slowness error matrix
    D=zeros(n_evt*n_sta,1);      % matrix of dt
    G=zeros(n_evt*n_sta,n_sta*3);   % make up matrix G(M*n_sta+1,n_sta+1)=0, operator matrix

    for i=1:n_evt
        for j=1:n_sta
            D( (i-1)*n_sta+j,    1    	) 	= sd(i,j) - sd(i,ref);
            
            G( (i-1)*n_sta+j,    3*j-2	)   = 1;
            G( (i-1)*n_sta+j,    3*j-1	)   = dx(i);
            G( (i-1)*n_sta+j,    3*j  	)   = dy(i);
        end
    end

%% calculate the Matrix M( n_sta*3,1  ) 
    M = G\D;      % Least Squares to solve  M
%     dSxy=M(1:3*n_sta);
    dS0=M(1:3:3*n_sta-2);
%     dS0=dS0-mean(dS0);   %added by liuwei, to eliminate the mean value
    dSx=M(2:3:3*n_sta-1);
    dSy=M(3:3:3*n_sta);

%% return dS
    dS2Dplus.lat_sta = lat_sta;
    dS2Dplus.lon_sta = lon_sta;
    dS2Dplus.dS0 = dS0;
    dS2Dplus.dSx = dSx;
    dS2Dplus.dSy = dSy;