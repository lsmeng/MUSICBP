% The Multiple Signal Classification (MUSIC) technique
%is a high-resolution technique designed
%to resolve closely spaced simultaneous sources. 
%MUSIC enables back-projection imaging (BP)
%with superior resolution than the standard Beamforming techinque.
%MUSICBP is a tutorial code in Matlab
%that images the spatial-temporal evolution of high-frequency radiators 
%(as a proxy of rupture front) for large earthquakes.
% You can modify simulation parameters in parts of the code 
% that have comments starting as: "**** Set here ..."
% September 3, 2019
% Lingsen Meng (meng@epss.ucla.edu) and Han Bao (hbrandon@ucla.edu)
clear
close all;
[status,address] =system('pwd'); addpath(address);   % add current directory to the path
%% *** Set here the processing steps to perform (positive=1, negtive=0)****
Initial_flag=0;     % Initializing a new project
readBP_flag=0;      % Reading seismogram from .SAC files
alignBP_flag=0;     % Hypocenter alignment
runBPbmfm_flag=0;   % Beamforming Back-projection
runBPmusic_flag=1;  % MUSIC Back-Projection

%% *** Set here the parameters to initialize the project and read the SAC files***
project = 'Palu_2018';% name of the project, e.g. Tohoku_2011
lon0=119.840;      	% hypocenter longitude
lat0=-0.178;      	% hypocenter latitude
dep=10.0;       	% hypocenter depth
Mw=7.5;            	% magnitude
sr=10;            	% sampling rate in Hz (the frequency that seismograms are down-sampled to) 
ori=60;             % length of seismograms before P-arrival time in seconds 
displayLength=360;  % length of waveforms (in seconds) to be displayed  
plotScale=1.5;      % amplitude scaling factor of seismograms for display purpose

%% *** Set here the Parameters for Hypocenter alignment ***
bandChoice=4;       % Choice of the alignment frequency band.
align(1,:)=[54,40,0.7];	% 1st align: freq band=[0.1, 0.25](Hz) windowLength=30(sec) maxShift=5(sec)
align(2,:)=[61,40,0.6]; % 2nd align: freq band=[0.25,0.5]  windowLength=15 maxShift=0.6
align(3,:)=[64,40,0.6]; % 3rd align: freq band=[0.5, 1.0]  windowLength=8  maxShift=0.1
align(4,:)=[64,0, 0.6];	% 4th align: freq band=[0.5, 1.0]  windowLength=8  maxShift=0.1
ts  = align(bandChoice,1); % start of the alignment window
refSta = align(bandChoice,2); % No. of the reference seismogram, set to zero for the stacked seismogram
cutoff= align(bandChoice,3); % cutoff threshold of the cross-correlation coefficient


%% *** Set here the Parameters for back-projection runner ***
inputBand=4;        % number of aligned seismograms used as an input for the back-projection
parr=61;            % P-wave arrival time, could be the same as ts
duration=45;        % earthquake
qs=60;              % number of grids in latitude of the imaging domain
ps=60;              % number of grids in longitude of the imaging domain
latrange=[-2 0.5];  % latitude length of the imaging domain
lonrange=[-0.5 0.5];% longitude length of the imaging domain
Band=4;             % frequency band for the back-projection
% Band=1 [0.05,0.25](Hz); Band=2 [0.25,1.0]; Band=3 [0.5,1]; Band=4 [0.5,2]; Band=5 [1,4]

%% *** NO CHANGE BELOW THIS LINE ***
workPath = './';
address = './';
path = [workPath,project,'/'];
if Initial_flag==1
    fprintf( ' Initializing... \n')
    initialize_BP(address,workPath,strcat(project,'/'));
end
if readBP_flag==1
	fprintf('Reading seismograms\n');
    Preshift=false; % no change!
    addfilelist(path)
    readteleBP(path,lon0,lat0,sr,ori,displayLength,Preshift,plotScale)
end
if alignBP_flag==1
    fprintf( 'Aligning Seismograms\n');
    check_data(path,1);
    align_BP(path,bandChoice,ts,refSta,cutoff,plotScale)
end
if runBPbmfm_flag==1
    fprintf( 'Running Beamforming back-projection\n');
    [BPfile]=setparBP(path,inputBand,Band,lon0,lat0,dep,parr,duration,ps,qs,lonrange,latrange);
    runteleBPbmfm(path,BPfile);
    summaryBPbeamforming
end
if runBPmusic_flag==1
    fprintf( 'Running MUSIC back-projection\n');
    [BPfile]=setparBP(path,inputBand,Band,lon0,lat0,dep,parr,duration,ps,qs,lonrange,latrange);
    runteleBPmusic(path,BPfile);
    summaryBPmusic
end

