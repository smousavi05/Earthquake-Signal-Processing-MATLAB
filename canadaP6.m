
clear all 
close all

%% Set the time range of interest
startTime = datenum(1989,12,29); % start time using MATLAB's datenum format
endTime = datenum(2017,10,30); % end time, 2 days after start time

%% Set the stations/channels to load
scnl = scnlobject('*', '*', '*', '*');

% Parameters for Calculate the wavelet transform -
opt.type = 'bump';           % Mother wavelet type 'bump'; 
opt.padtype = 'symmetric';   % padded via symmetrization
opt.rpadded = 1;
opt.nv = 8;                 % Number of voices

datapath = '../../../../../../Volumes/MyBook/canadaTomo/canadaNew5';
outpath = '../../../../../../Volumes/MyBook/canadaTomo/canadaNew6';
all_events = dir(datapath);
events = all_events(4:length(all_events));
num_dir = numel(events);
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
for evnum=1:numel(events); 
% for evnum=10:12;
events(evnum).name
sacfiles = dir(fullfile(datapath, events(evnum).name, '*.SAC'));

%% Set the datasource for each file
for filenum=1:5
%for filenum=1:length(sacfiles)
    ds = datasource('sac', fullfile(datapath, events(evnum).name, sacfiles(filenum).name) );
    w = waveform(ds, scnl, startTime, endTime);
    
  arclen = distance(get(w,'EVLA'),get(w,'EVLO'),get(w,'STLA'),get(w,'STLO'),'degrees');
  dss = arclen.*111120
%     dis = get(w,'DIST')
%     dd = dis.*1000;
 w = addfield(w,'distance', dss);
 w = addfield(w,'KNETWK', 'CN');    
savesac(w,fullfile(outpath, events(evnum).name,sacfiles(filenum).name));


end    
end

