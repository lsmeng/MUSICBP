% USAGE: after download and rdseed all events in mother_loc, this script
%        will:
%           1. automatically read all .SAC files for each event.
%           2. match stations with MAINSHOCK.
%           3. transport 'timeshift' to matched stations.
%           4. Generate data_matched.mat with matched stations and
%           timeshift(mainshock alignment)
% INPUT: 1. mother_loc: Mother directory of all events
%        2. main_mat: .mat file of the Mainshock
%        3. ori:      for aftershock, how many sec before P-arrival
%        4. display:  for aftershock, how many sec of signal in total (before_P + after_P)

clear
close all

%% INPUT
mother_loc = '/home/liuwei/HPSSD/Qinghai/Tele_Cali/TeleDD_Loc/AU/Events';
main_mat   = '/home/liuwei/HPSSD/Qinghai/BP/TeleDD_Loc/AU/Input/Par0.5_2_10.mat';
address= '/home/liuwei/HPSSD/Supershear/BP2/Workshop2019/';            % address of paths that need to be added below
addpath(genpath(strcat(address,'funcLib/')));   % directory of Back-Projection functions
ori      = 60; % please refer to your request on IRIS
display  = 360; % please refer to your request on IRIS
Preshift = false;	% no change! About readBP
plotscale= 1.5;     % plotscale:  scaling of the amplitudes of seismograms


%% Loop for each aftershock
cd (mother_loc)
load evtlst.mat
[n_evt,~] = size(evtlst_name(:,1)); % n_evt: number of aftershock events

for i = 1 : n_evt 
   %%% Load mainshock .mat file
    load(main_mat);
    ret_main = ret;
   
   %%% readBP for aftershock events
    path = [mother_loc '/' deblank(evtlst_name(i,:)) '/'];
    cd (path)
    mkdir ('Input')
    cd ('Data')
    ! ls *.SAC > filelist
    cd ..
    readteleBP(path,ret_main.lon0,ret_main.lat0,ret_main.sr,ori,display,Preshift,plotscale)
    
    %%%check for data
    %check_data(path,1);
    
    %%% Match stations with mainshock 
    load('Input/data0.mat')
    [~,i_main,i_aft] = intersect(deblank(string(ret_main.nm)), strtrim(string(ret.nm))); % find same station
        %%%% Since the 'recordtime' of different event downloaded from IRIS
        %%%% might be different, we also need to correct this difference,
        %%%% by transit this error into aftershock's 'timeshift'
        ret = orderAll(ret,0,i_aft);
        ret_main = orderAll(ret_main,0,i_main);
        rctime_diff_sec      = 3600*24*(ret.recordtime      - ret.recordtime(1));
        rctime_main_diff_sec = 3600*24*(ret_main.recordtime - ret_main.recordtime(1));
        diff = rctime_main_diff_sec - rctime_diff_sec;
        ret.timeshift = ret_main.timeshift + diff;
        for j = 1:length(diff)
            ret.xori(j,:) = specshift(ret.xori(j,:),ret.timeshift(j)*ret.sr);
        end
    
    save Input/data5.mat ret
end
%%% TEST


% tmp1 = (ret.recordtime*3600*24 + ret.timeshift) - (ret.recordtime(1)*3600*24 + ret.timeshift(1));
% tmp2 = (ret_main.recordtime*3600*24 + ret_main.timeshift) - ...
%     (ret_main.recordtime(1)*3600*24 + ret_main.timeshift(1));
% tmp1-tmp2
% x1(i,:)=specshift(x0(i,:),timeshift(i)*sr);

    