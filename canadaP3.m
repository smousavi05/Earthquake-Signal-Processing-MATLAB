clear all
close all

%% defining the parameters
minla = 40;  % minimumlatitude
maxla = 65;  % maximumlatitude
minlo = -85; % minimumlongitude
maxlo = -40  % maximumlongitude

minmag = 1.5; % minimummagnitude

starttime = '1989-01-21';
endtime = '2017-01-20';


%% finding events
%run GISMO/startup_GISMO.m;
ds = datasource('irisdmcws');

ev = irisFetch.Events('boxcoordinates', [minla, maxla, minlo, maxlo] ...
, 'starttime', starttime ,'endtime',endtime,'minimummagnitude',minmag);

% cobj = Catalog.retrieve('iris', 'boxcoordinates', [minla, maxla, minlo, maxlo] ...
% , 'starttime', starttime ,'endtime',endtime,'minimummagnitude',minmag);

fileID = fopen('CanadaCat.txt','w');

for i=1:length(ev)
fprintf(fileID,'%2.4f %2.4f %2.4f %1.2f %6s\n',ev(i).PreferredLatitude,ev(i).PreferredLongitude ...
    ,ev(i).PreferredDepth,ev(i).PreferredMagnitudeValue,ev(i).PreferredMagnitudeType);
end
fclose(fileID);

% evtr =[];
% ds = datasource('irisdmcws');
% for ke = 1:length(cobj.otime);
%     disp(sprintf('Event Number: %d',ke))
%      et = cobj.otime(ke);t1 = et - minutes(5); t2 = et + minutes(12);
%      sTime = datestr(t1,'yyyy-mm-dd hh:MM:ss'); eTime = datestr(t2,'yyyy-mm-dd hh:MM:ss');
%      
%      % finding available network at the time of each event
%      arcticNetworkList = irisFetch.Networks('Station','','','','','boxcoordinates' ...
%      , [minla, maxla, minlo, maxlo],'starttime', sTime,'endtime',eTime);
%      
%        for kn = 1:length(arcticNetworkList);
%             for ks = 1:length(arcticNetworkList(kn).Stations);
%                 
%                  this_sta = char(arcticNetworkList(kn).Stations(ks).StationCode);
%                  this_net = char(arcticNetworkList(kn).NetworkCode);
%                  this_chan = '?HZ';
%                  
%                  % Here you can use scnlobject or the newer ChannelTag -
%                  % waveform understands both
%                  scnl = scnlobject( this_sta, this_chan, this_net, '--');
%                  %scnl = ChannelTag(this_net, this_sta, '--', this_chan);
%                 
%                 w = waveform(ds, scnl, sTime, eTime);
%                 if ~isempty(w)
%                     plot_panels(w, true);
%                     evtr = [evtr w];
%                 end
%           end
%        end
%        
%     disp('Number of Traces:')
%     disp(length(evtr))
% end 
