
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

datapath = '../../../../../../Volumes/MyBook/canadaTomo/canadaNew6';
outpath = '../../../../../../Volumes/MyBook/canadaTomo/canadaNew5';
all_events = dir(datapath);
events = all_events(4:length(all_events));
num_dir = numel(events);
%%%%%%%%%%%%%%%%%%%%%%
for evnum=1:numel(events); 
% for evnum=10:12;
events(evnum).name
sacfiles = dir(fullfile(datapath, events(evnum).name, '*.SAC'));

%% Set the datasource for each file 
for filenum=1:length(sacfiles)
    ds = datasource('sac', fullfile(datapath, events(evnum).name, sacfiles(filenum).name) );
    w = waveform(ds, scnl, startTime, endTime);

    fr = get(w,'freq'); dis = get(w,'DIST');
if fr >= 20 & dis >= 150
%plot_panels(w, true);

%% Clean data
w = addfield(w,'Parriv', get(w,'T1'));
w = addfield(w,'DELTA', 1/get(w,'freq'));
w = medfilt1(w,3); % remove any spikes of sample length 1
w = fillgaps(w, 'interp');
w = detrend(w);
w = demean(w); %remove the voltage offset from the data.

%% Specify polezero file
respath = '../../../../../../Volumes/MyBook/canadaTomo/responceCanada/';
v = sprintf('%s%s%s%s%s%s%s', respath, 'SAC_PZs_',get(w,'network'),'_', get(w,'station'), '__', get(w,'channel'));
pz = sacpz(sprintf(v));

% tw =datestr(get(w,'start'), 'yyyy mm dd HH MM')

if length(pz) == 1
    
  polezero.poles = pz.p;
  polezero.zeros = pz.z;
  polezero.normalization = pz.k;

else
pzt=zeros(length(pz),2);
for i=1:length(pz)
    pzt(i,1) = pz(i).starttime;
    pzt(i,2) = pz(i).endtime;
end

[c] = find(pzt(:,1)<= get(w,'start')& pzt(:,2)>= get(w,'start'))
if length(c) == 0
  polezero.poles = pz(end).p;
  polezero.zeros = pz(end).z;
  polezero.normalization = pz(end).k
else
polezero.poles = pz(c).p;
polezero.zeros = pz(c).z;
polezero.normalization = pz(c).k;
end
end

% frequencies = logspace(-2,2,100);
% response_polezero = response_get_from_polezero(frequencies,polezero);
% response_plot(response_polezero)

%% Deconvolve instrument response 
filterObj = filterobject('b',[0.4 9],3);
wFilt = filtfilt(filterObj,w);
wCorrected = response_apply(wFilt,filterObj,'polezero',polezero);
%plot_panels(wCorrected)

snan = sum(isnan(get(wCorrected,'data')));
 if snan < 10;
ss = strsplit(events(evnum).name,'_');
OO = strsplit(events(evnum).name,'.');
eVname =  sprintf('%s.%s.%s.M%s',char(OO(1)),char(OO(2)),char(OO(3)),char(ss(2)))
SP = sprintf('%s%s', '../../../../../../Volumes/MyBook/canadaTomo/canadaNew5/',eVname);

sTname1 = sprintf('IC.%s.%s.%s.%s.SAC', eVname, get(w,'network') ...
              , get(w,'station'),get(w,'channel'));
sTname2 = sprintf('DN.%s.%s.%s.%s.SAC', eVname, get(w,'network') ...
              , get(w,'station'),get(w,'channel'));
          
 %% denoising
    noisy1 = get(wCorrected,'data');
    [wnoisy1,as] = cwt_fw(noisy1,opt.type,opt.nv,get(wCorrected,'DELTA'));
    [na n] = size(wnoisy1);
    at = get(w,'freq')* get(w,'T1');
                   
    % Round 1
      for L = 1:na
          if length(wnoisy1) < at; at = length(wnoisy1)/3; end
                   
          if any(wnoisy1(L,:)) == 1;
             T = rms(abs(wnoisy1(L,1:at)));
           for R = 1:n;   
              if abs(wnoisy1(L,R)) <= T
                 wnoisy1(L,R) = 0;
              else
                 res = abs(wnoisy1(L,R)) - T;
                 res = (res + abs(res))/2;
                 wnoisy1(L,R) = sign(wnoisy1(L,R))*res;
               end
             end
          end   
        end
        dn1 = cwt_iw(wnoisy1, opt.type, opt);
                 
    % Round2 
    [wnoisy2,as2] = cwt_fw(dn1,'morlet',opt.nv,get(wCorrected,'DELTA'));
      for k = 1:na
         if  any(wnoisy2(k,:)) == 1;
             Wx_fine = abs(wnoisy2(k, 1:at));
             lamba = sqrt(2*log(at)) * mad( abs(Wx_fine (:))) * 1.4826;
             wnoisy2(k,:) = wnoisy2(k,:).* (abs(wnoisy2(k,:)) > lamba);
         end
      end
      dn2 = cwt_iw(wnoisy2, opt.type, opt);

%% saving the results
A = exist(sprintf(SP));

if A == 0
     mkdir('../../../../../../Volumes/MyBook/canadaTomo/canadaNew5/',eVname);
end

% saving instrument corrected signal
savesac(wCorrected,fullfile(outpath, eVname),sTname1);

% saving denoised signal
wCorrected = set(wCorrected,'data', dn2); % was traces(i).getData();
wCorrected = medfilt1(wCorrected,10); % remove any spikes of sample length 1
savesac(wCorrected,fullfile(outpath, eVname),sTname2);
%plot_panels(wCorrected, true)                 
% 
% %% Plot result
% plot_panels(wCorrected)
% % plot_spectrum(w);
% % spectrogram(w);
% 
close all
end    
end
end
end