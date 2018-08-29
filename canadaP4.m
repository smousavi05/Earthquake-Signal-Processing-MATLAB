clear all
close all

%% defining the parameters
minla = 40;  % minimumlatitude
maxla = 65;  % maximumlatitude
minlo = -85; % minimumlongitude
maxlo = -40; % maximumlongitude

minmag = 2.5; % minimummagnitude

maxDepth = 35; % maximum depth of events

starttime = '2000-01-01';
endtime = '2004-02-02';

mindist = 1.5; % minimum distance between station and event (degree)
maxdist = 30; % maximum distance between station and event (degree)

minfrq = 40; % minimum sample rate of records

X1 = ['Minimum Latitude:',num2str(minla),'  Maximum Latitude:',num2str(maxla)]; disp(X1)
X2 = ['Minimum Longitude:',num2str(minlo),'  Maximum Longitude:',num2str(maxlo)]; disp(X2)
X3 = ['Start Time:',starttime,'  End Time:',endtime];disp(X3)
X4 = ['Minimum Distance:',num2str(mindist),'  Maximum Distance:',num2str(maxdist)]; disp(X4)
X5 = ['Minimum Magnitude:',num2str(minmag)];disp(X5)
X6 = ['Maximum Depth:',num2str(maxDepth)];disp(X6)

% Parameters for Calculate the wavelet transform -
opt.type = 'bump';           % Mother wavelet type 'bump'; 
opt.padtype = 'symmetric';   % padded via symmetrization
opt.rpadded = 1;
opt.nv = 8;                 % Number of voices

%% finding events
run GISMO/startup_GISMO.m;
ds = datasource('irisdmcws');

ev = irisFetch.Events('boxcoordinates', [minla, maxla, minlo, maxlo] ...
, 'starttime', starttime ,'endtime',endtime,'minimummagnitude',minmag,'maximumDepth',maxDepth);

eV =[]; sV =[]; trnum = 0;
for ke = 1:numel(ev);
   X7 = ['Event Number: ', num2str(ke),' from ', num2str(numel(ev)), ' events']; disp(X7)
   
    
     et = datenum(ev(ke).PreferredTime);t1 = et - minutes(5); t2 = et + minutes(12);
     sTime = datestr(t1,'yyyy-mm-dd hh:MM:ss'); eTime = datestr(t2,'yyyy-mm-dd hh:MM:ss');
     
     % finding available network at the time of each event
     arcticNetworkList = irisFetch.Networks('Station','','','','','boxcoordinates' ...
     , [minla, maxla, minlo, maxlo],'starttime', sTime,'endtime',eTime);
     
       for kn = 1:numel(arcticNetworkList);
            for ks = 1:numel(arcticNetworkList(kn).Stations);
                
                         
            traces = irisFetch.Traces(char(arcticNetworkList(kn).NetworkCode) ...
                ,char(arcticNetworkList(kn).Stations(ks).StationCode) ...
                ,'*','?HZ',sTime, eTime,'includePZ');
            
          
            if isempty(traces) == 0;
                  for kt = 1:length(traces)
                   close all
                  trnump = trnum + 1; trnum = trnump; 
                  X8 = ['Total Number of Traces: ', num2str(trnum)]; disp(X8)
   
                      
                  [arclen,az] = distance(ev(ke).PreferredLatitude,ev(ke).PreferredLongitude ...
                      ,traces(kt).latitude,traces(kt).longitude);
                  
                  [brclen,baz] = distance(traces(kt).latitude,traces(kt).longitude, ...
                      ev(ke).PreferredLatitude,ev(ke).PreferredLongitude);
                  
                  if traces(kt).sampleRate >= minfrq & mindist <= arclen <= maxdist ;
                  clear w
                  w = waveform;
                  chaninfo = ChannelTag(char(traces(kt).network), ...
                  char(traces(kt).station), ...
                  char(traces(kt).location), ...
                  char(traces(kt).channel));
                  w = set(w,'channelinfo',chaninfo,'freq',traces(kt).sampleRate); %, 'start', datenum(startDateStr, 'yyyy-mm-dd HH:MM:SS.FFF'));
                  w = set(w,'start', char(datestr(traces(kt).startTime, 'yyyy-mm-dd HH:MM:SS')));
                  w = addfield(w,'DELTA', 1/traces(kt).sampleRate);
                  w = addfield(w,'STLA', traces(kt).latitude);
                  w = addfield(w,'STLO', traces(kt).longitude);
                  w = addfield(w,'STEL', traces(kt).elevation);
                  w = addfield(w,'STDP', traces(kt).depth);
                  w = addfield(w,'EVLA', ev(ke).PreferredLatitude);
                  w = addfield(w,'EVLO', ev(ke).PreferredLongitude);
                  w = addfield(w,'EVDP', ev(ke).PreferredDepth);                 
                  w = addfield(w,'EVMAG', ev(ke).PreferredMagnitudeValue);
                  w = addfield(w,'AZ', az);
                  w = addfield(w,'BAZ', baz);
                  w = addfield(w,'dip', traces(kt).dip);
                  w = addfield(w,'sensitivity',traces(kt).sensitivity);
                  w = addfield(w,'sensitivityFrequency',traces(kt).sensitivityFrequency);
                  w = addfield(w,'instrument',char(traces(kt).instrument));
                  w = set(w,'units',char(traces(kt).sensitivityUnits));
                  w = addfield(w,'calib',1 ./ traces(kt).sensitivity);
                  w = addfield(w,'calib_applied','NO');
                  w = set(w,'data', traces(kt).data); % was traces(i).getData();
                 
                  if get(w,'duration') > minutes(2)
                  
                  %% Clean data
                  w = medfilt1(w,3); % remove any spikes of sample length 1
                  w = fillgaps(w, 'meanAll');
                  w = detrend(w);
                  w = demean(w); %remove the mean.
     
                 %% Deconvolve instrument response 
                  clear polezero.poles; clear polezero.zeros; clear polezero.normalization;
                  polezero.poles = traces(kt).sacpz.poles;
                  polezero.zeros = traces(kt).sacpz.zeros;
                  polezero.normalization = traces(kt).sacpz.constant;
 
                  filterObj = filterobject('b',[0.5 15],2);
                  wFilt = filtfilt(filterObj,w);
                  wCorrected = response_apply(wFilt,filterObj,'polezero',polezero);
                  
                  %% check to see if any event exist in the trace
                  Z = sta_lta(wCorrected);  
                  
                  
                  if Z.evn >= 1 & Z.evn < 5 ;
                      
                  startTrace = datestr(get(wCorrected, 'start'), 'mmmm dd, yyyy HH:MM:SS.fff');
                  t1 = datevec(startTrace,'mmmm dd, yyyy HH:MM:SS.fff');
                  
                  startT = datestr(Z.stTrig(1,:), 'mmmm dd, yyyy HH:MM:SS.fff');
                  t2 = datevec(startT,'mmmm dd, yyyy HH:MM:SS.fff');
                  trigTime = etime(t2,t1)
                  
                  stDN = datenum(Z.stTrig(1,:));
                  correctedT = stDN - minutes(2); tC = datestr(correctedT,'yyyy-mm-dd hh:MM:ss');
                  triPlus =  stDN + minutes(5); tP = datestr(triPlus,'yyyy-mm-dd hh:MM:ss');
                  triMinus = datenum(tC) - minutes(5); tM = datestr(triMinus,'yyyy-mm-dd hh:MM:ss');
                  wP = extract(wCorrected, 'TIME', startT, tP); RMSP = rms(wP);                 
                  wM = extract(wCorrected, 'TIME', tM, startT); RMSM = rms(wM);
                  snr = RMSP/RMSM;
                  
                  if snr >= 2                
                
                  close all                  
                  Ptime = PphasePicker(get(wCorrected,'data'), get(wCorrected,'DELTA'), 'wm', 'Y', 0.01, 0.7, 400, 'to_peak')
                   
                   if abs(Ptime - trigTime) < 50;
                   w = addfield(w,'Parriv', Ptime);
                   else
                   w = addfield(w,'Parriv', trigTime);
                   end
                   get(w,'Parriv')
                   
                   %% denoising
                   noisy1 = get(wCorrected,'data');
                   [wnoisy1,as] = cwt_fw(noisy1,opt.type,opt.nv,get(wCorrected,'DELTA'));
                   [na n] = size(wnoisy1);
                   at = traces(kt).sampleRate * get(w,'Parriv');
                   
                   % Round 1
                   for L = 1:na
                       if  any(wnoisy1(L,:)) == 1;
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
                                      
                   eVname = sprintf('%s.M%s', datestr(ev(ke).PreferredTime,'yyyy-mm-ddThh.MM.ss') ...
                      , num2str(ev(ke).PreferredMagnitudeValue))
                   outpath = ('../../../../../../Volumes/MyBook/canadaTomo/canadaNew4');
                   mkdir(outpath,eVname);
                 
                   % saving instrument corrected signal
                   sTname = sprintf('IC.%s.%s.%s.%s.SAC', eVname, get(w,'network') ...
                      , get(w,'station'),get(w,'channel'))
    
                   savesac(wCorrected,fullfile(outpath, eVname),sTname);
                   
                   % saving denoised signal
                   wCorrected = set(wCorrected,'data', dn2); % was traces(i).getData();
                   wCorrected = medfilt1(wCorrected,10); % remove any spikes of sample length 1
                   
                   sTname = sprintf('DN.%s.%s.%s.%s.SAC', eVname, get(w,'network') ...
                      , get(w,'station'),get(w,'channel'))
    
                   savesac(wCorrected,fullfile(outpath, eVname),sTname);
                   plot_panels(wCorrected, true)              
%                    pause
                   e = sprintf('%s %7s %8s %4s %4s %4s', datestr(ev(ke).PreferredTime,'yyyy-mm-ddThh.MM.ss') ...
                        , num2str(ev(ke).PreferredLatitude) ...
                        , num2str(ev(ke).PreferredLongitude) ...
                        , num2str(ev(ke).PreferredDepth) ...
                        , num2str(ev(ke).PreferredMagnitudeValue) ...
                      , num2str(ev(ke).PreferredMagnitudeType)); 
               
%                       eV =[eV; e]
                      
                  es = sprintf('%4s %5s %8s %8s', get(w,'network') ...
                        , get(w,'station') ...
                        , num2str(get(w,'STLA')) ...
                        , num2str(get(w,'STLO'))...
                        ); 
                    
%                       sV =[sV; es]
                  end
                   end
                 end
             end
          end   
      end
  end
  end
end

% fileID = fopen('cat.txt','w');
% C = unique(eV,'rows')
% fprintf(fileID,'%s %s %s %s %s %s \n',C');
% fclose(fileID);

