% USAGE: after bb201, we've already had the raw data for each event along
%        with the corresponding timeshift of each station given by mainshock.
%        Now, this script will assign parameters for runBP, and calculate
%        beamforming backprojection in such a period that containing each 
%        aftershock event (30s before P-arrival to 30s after P-arrival).
%        This script will finally generate a 'logfile' for each event that
%        consists with 4 columns, e.g. 
%               10.000000 37.708696 145.860000 10287.839517
%        which are relative time, lat, lon, and relative power.
% 
% INPUT: 1. mother_loc: Mother directory of all events
%        2. evtlst: event list
%        3. .mat data for runBP (seismic data)


clear
close all
%% INPUT
mother_loc = '/home/liuwei/HPSSD/Qinghai/Tele_Cali/TeleDD_Loc/AU/Events';
address='/home/liuwei/HPSSD/Supershear/BP2/Workshop2019/';            % address of paths that need to be added below
addpath(genpath(strcat(address,'funcLib/')));   % directory of Back-Projection functions	

%%
cd (mother_loc)
load 'evtlst.mat'; % load event list
[n_evt,~] = size(evtlst_char);

%% Loop for all events
% for j=1:2
for j = 1 : n_evt
    
    load 'evtlst.mat'
    evtlst_name(j,:)
    cd (evtlst_name(j,:))
        mother_loc = '/home/liuwei/HPSSD/Qinghai/Tele_Cali/TeleDD_Loc/AU/Events';
        % Parameter set for runBP
            lon0=98.3622588499; 	% MAINSHOCK lontitude
            lat0=34.6202641954;  	% MAINSHOCK latitude
            dep0 = 7.607;    % MAINSHOCK depth
            sr=10;          % sample rate
            parr=45;      	% start time
            begin=0;       	% always 0
            over=30;    	% end time
            step=1;         % time step
            ps=120;          % number of grids for lat
            qs=120;          % number of grids for lon
            latrange=[-0.8,0.8]; % lat range
            lonrange=[-1.5,1.5]; % lon range
            fl=0.5;        % frequence low of bandpass
            fh=2;           % frequence high of bandpass
            win=8;         % window length for BP

            inputband=5;    % i.e. data5.mat
            Band = 4; % band4: [0.5,2]
    path=strcat(mother_loc,'/',evtlst_name(j,:),'/');    

    [BPfile]=setparBP(path,inputband,Band,lon0,lat0,dep0,parr,over,ps,qs,lonrange,latrange);
    runteleBPmusic(path,BPfile)
    posteriorBPmusic;
    close all;
    % END of current loop
    cd ../../..
end
fclose(fid1);
