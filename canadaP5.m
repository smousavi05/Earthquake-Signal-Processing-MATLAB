clear all
close all
%% Checking Code
%% finding events
all_events = dir('../../../../../../Volumes/MyBook/canadaTomo/canadaNew5/*');
events = all_events(4:length(all_events));
num_dir = numel(events)
prompt = 'Delete(1)? ';


%%%%%   for manual quality checking 
 for ii=1:num_dir
     
    disp(events(ii).name)
    v = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5', num2str(events(ii).name),'DN.*');
    all_seismograms = dir(sprintf(v));
    num_seis = numel(all_seismograms);
    
    vv = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5', num2str(events(ii).name),'IC.*');
    all_seismograms2 = dir(sprintf(vv));

    
    for jj = 1:num_seis
        disp(all_seismograms(jj).name)
        comName = strsplit(sprintf('%s',all_seismograms(jj).name),'.');
        if length(comName) == 10
        seisName = strcat('*.',comName(7),'.',comName(8),'.',comName(9),'.','*');
        snameIC =  strcat('IC.*.',comName(7),'.',comName(8),'.',comName(9),'.','*');
        elseif length(comName) == 9
        seisName = strcat('*.',comName(7),'.',comName(8),'.','*');
        end
        
        P = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5', num2str(events(ii).name),char(seisName));
        B = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5', num2str(events(ii).name),num2str(all_seismograms(jj).name));
        BB = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5', num2str(events(ii).name),num2str(all_seismograms2(jj).name));
        [tD,dataD,hdr] = read_sac(sprintf(B));
        [tI,dataI,hdr] = read_sac(sprintf(BB));
%%       
        figure 
        subplot 211
         plot(tD,dataD)
        subplot 212
         plot(tI,dataI)
         
         str = input(prompt);
         if str == 1; delete(sprintf(P));end

% %
% dt = hdr.times.delta;
% t = linspace(0,(hdr.times.e - hdr.times.b),length(data));
% [a indx]=find(data == max(data));
% ts = a(1)*dt - 10;
% te = a(1)*dt + 100;
% tn = ts - 110; 
% if te/dt > length(data); te = length(data)*dt; end
% if tn < 1; tn = 1;end
% if ts < 0; delete(sprintf(P));
% else
% RMSS = rms(data(ts/dt:te/dt));
% RMSN = rms(data(tn/dt:ts/dt));
% snr = RMSS/RMSN
% if snr < 2; delete(sprintf(P));end
% end

% [arclen,az] = distance(hdr.event.evla,hdr.event.evlo, ...
%                   hdr.station.stla,hdr.station.stlo);
%                   
%  if arclen < 1.5; delete(sprintf(P));end  
        
       close all
    end
 end

 


% %%  
%   for ii=1:num_dir
%       v = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5', num2str(events(ii).name),'DN.*');
%       all_seis = dir(sprintf(v));
%       num = numel(all_seis);
%      if  num < 3
%          disp(events(ii).name)
%         rmdir(sprintf('%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5', num2str(events(ii).name)),'s');
%      end
%   end
  
% %%  for making events.par file as Zhao's format
% fid=fopen('events.txt','w');
% fprintf(fid,'dir                        year mm dd hh min s Lat.     Lon.        Dep.      m\n');
% 
% for ii=1:num_dir
%     v = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5',...
%         num2str(events(ii).name),'DN.*');
%     all_seismograms = dir(sprintf(v));
%     B = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5',...
%         num2str(events(ii).name),num2str(all_seismograms(1).name));
%     [t,data,hdr] = read_sac(sprintf(B));
% 
% evName = strsplit(sprintf('%s',events(ii).name),'-');
% evName2 = strsplit(sprintf('%s',events(ii).name),'T');
% evName3 = strsplit(sprintf('%s',char(evName2(1))),'-');
% evName4 = strsplit(sprintf('%s',char(evName2(2))),'.');
% evName5 = strsplit(sprintf('%s',events(ii).name),'M');
% 
%  fprintf(fid,'%25s  %4s %2s %2s %2s %2s %2s %4f %4f %2f  %3s \n',events(ii).name,char(evName(1)),...
%      char(evName(2)),char(evName3(3)),char(evName4(1)),char(evName4(2))...
%      ,char(evName4(3)),hdr.event.evla,hdr.event.evlo,hdr.event.evdp,char(evName5(2)));
% 
% end
% fclose(fid);
 


% %%  for making stats.d file for Zhao's format 
% fid=fopen('Oklahoma2.stats.txt','w');
% fprintf(fid,'Lat.      Lon.      Sta.     Network\n');
% for ii=1:num_dir
%     v = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5',...
%         num2str(events(ii).name),'DN.*');
%     all_seismograms = dir(sprintf(v));
%     num_seis = numel(all_seismograms);
%     for jj = 1:num_seis
% 
%     B = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5',...
%         num2str(events(ii).name),num2str(all_seismograms(jj).name));
%     [t,data,hdr] = read_sac(sprintf(B));
% 
% stName = strsplit(sprintf('%s',all_seismograms(jj).name),'.');
% lo = num2str(hdr.station.stlo);
% la = num2str(hdr.station.stla)
% 
% if length(stName) == 10; st = char(stName(8)); nt = char(stName(7));end
% if length(stName) == 9; st = char(stName(7)); nt = char(stName(6));end
% 
%  fprintf(fid,'%4s  %4s   %5s    %2s \n',lo,la,st,nt);
%     end
%  
% 
% end
% fclose(fid);
 
% % cat Canada.stats.txt |  sort -u > stat.txt
%




% %% for making event list for GMT plotting
% fid=fopen('Oklahoma2.events.txt','w');
% for ii=1:num_dir
%     v = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5',...
%         num2str(events(ii).name),'DN.*');
%     all_seismograms = dir(sprintf(v));
%     B = sprintf('%s/%s/%s','../../../../../../Volumes/MyBook/canadaTomo/canadaNew5',...
%         num2str(events(ii).name),num2str(all_seismograms(1).name));
%     [t,data,hdr] = read_sac(sprintf(B));
% 
% evName = strsplit(sprintf('%s',events(ii).name),'-');
% evName2 = strsplit(sprintf('%s',events(ii).name),'T');
% evName3 = strsplit(sprintf('%s',char(evName2(1))),'-');
% evName4 = strsplit(sprintf('%s',char(evName2(2))),'.');
% evName5 = strsplit(sprintf('%s',events(ii).name),'M');
% 
%  fprintf(fid,'%4f %4f %2f  %3s \n',hdr.event.evlo,hdr.event.evla,hdr.event.evdp,char(evName5(2)));
% 
% end
% fclose(fid);
