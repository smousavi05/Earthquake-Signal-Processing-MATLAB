
clear all
close all

% % 1 : long
% % 2 : lat
% % 3 : year
% % 4 : month
% % 5 : day
% % 6 : mag
% % 7 : depth
% % 
% % % E
% % A = xlsread('from-site-declust','A');
% % [long lat yr mon day mag dpt h m ] = textread ('resenburg-whole.txt', '%f %f %f %f %f %f %f %f %f');
% % A = [long lat yr mon day mag dpt h m];
% % 
% % [r c]= size(A);
% % 
% % out=[];
% %  for d=30:0.5:80
% %     for mm = 10:0.5:45
% %      oo = [d, mm, 0];
% %      out = [out; oo];
% %     end
% %  end
% % 
% % oo=[];
% % 
% % for i=1:r;
% % 
% % if ((A(i,3) >= -1250) && (A(i,3) < 1900) && (A(i,6) >= 5.0))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% %  
% % if ((A(i,3) >= 1900) && (A(i,3) < 1927) && (A(i,6) >= 5.0))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% % elseif ((A(i,3) >= 1927) && (A(i,3) < 1963) && (A(i,6) >= 4.5))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% %  
% % elseif ((A(i,3) >= 1963) && (A(i,3) < 1973) && (A(i,6) >= 4.8))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% % elseif ((A(i,3) >= 1973) && (A(i,3) < 1995) && (A(i,6) >= 4.5))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% % elseif  ((A(i,3) >= 1995) && (A(i,3) < 2001) && (A(i,6) >= 4.2))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% % elseif ((A(i,3) >= 2001) && (A(i,3) < 2007) && (A(i,6) >= 4.0))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% %     end
% % end
% % 
% % 
% % %% magnitude /time distribution
% % plot(oo(:,3),oo(:,1),'o')
% % grid on
% % xlabel('Time (years)','FontSize',18,'FontWeight','bold')
% % ylabel('Longitude (degree)','FontSize',18,'FontWeight','bold')
% % set(gca,'FontSize',15)
% % 
% % % xlim([1995,2001])
% % % title('4.8 - 5.4')
% % 
% %   mm =sum(out(:,3))/length(out);
% %  out(:,3) = out(:,3)-mm;
% %  mino=min(out(:,3));
% %  
% %  for d=30:0.5:80
% %     for mm = 10:0.5:45
% %      oo = [d, mm, mino];
% %      out = [out; oo];
% %     end
% %  end
% % 
% %  
% % delete E.txt
% % fileID = fopen('E.txt','w');
% % fprintf(fileID,'%f %f %f %f\n',out'); %E
% % fprintf(fileID,'%f %f %f \n',out'); %D
% % fclose(fileID);
% % movefile('E.txt','../eng3/.') 
% % system('./eng3/runEng.sh')
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%
% % 1 : long
% % 2 : lat
% % 3 : year
% % 4 : month
% % 5 : day
% % 6 : mag
% % 7 : depth

%% E
[name time h min sec lat long dpth dd kk mag mt hh] = textread ('eventsLLLg.txt', '%s %f %f %f %f %f %f %f %f %f %f %s %f');
A = [long lat yr mon day mag dpt h m];
[r c]= size(A);
AA =[];

for i=1:r;
% if ((A(i,3) >= -1250) && (A(i,3) < 1900) && (A(i,6) >= 4.0))
% b =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; 
% AA = [AA;b];
% 
% elseif ((A(i,3) >= 1900) && (A(i,3) < 1927) && (A(i,6) >= 4.0))
% b =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; 
% AA = [AA;b];
% 
% elseif ((A(i,3) >= 1927) && (A(i,3) < 1963) && (A(i,6) >= 4.0))
% b =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; 
% AA = [AA;b];
% 
% elseif ((A(i,3) >= 1963) && (A(i,3) < 1973) && (A(i,6) >= 4.8))
% b =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; 
% AA = [AA;b];
% 
% elseif ((A(i,3) >= 1973) && (A(i,3) < 1995) && (A(i,6) >= 4.3))
% b =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; 
% AA = [AA;b];

if ((A(i,3) >= 1995) && (A(i,3) < 2001) && (A(i,6) >= 4.2))
b =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; 
AA = [AA;b];

elseif ((A(i,3) >= 2001) && (A(i,3) < 2007) && (A(i,6) >= 4.0))
b =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; 
AA = [AA;b];

end
end


% grid info
dx = 0.5
dy = 0.5
delta = sqrt(dx^2+dy^2);
out=[];

% Loop over Longitudes
 for xx=30:dx:80
% Loop over Latitudes
    for yy = 10:dy:45

    l = sqrt(((AA(:,1) - xx)).^2 + ((AA(:,2) - yy)).^2);
    ll = l(find(l <= delta)); 
    [tf, index] = ismember(ll, l);
    gv = sum(AA(index,3));
%    gv = mean(AA(index,3));
    AA(index,:)=[];
        
    o = [xx, yy, gv];
    out = [out; o];
    end
 end

 out(isnan(out(:,3)),:)=[];
 out(find(out(:,3) == 0),:)=[];
 out(:,3)=log(out(:,3));
 
 
delete E.txt
fileID = fopen('E.txt','w');
fprintf(fileID,'%.2f %.2f %.2f \n',out'); %E
fclose(fileID);



% dx = 0.5;
% dy = 0.5;
% 
% minlon=30;
% maxlon=80;
% minlat=10;
% maxlat=45;
% nx =(maxlon-minlon)/dx;
% ny =(maxlat-minlat)/dy;
% 
% xg = linspace( minlon, maxlon, nx+1 ) ;
% yg = linspace( minlat, maxlat, ny+1 ) ;
% nCells = nx * ny ;
% 
% figure(1) ;  clf ;  hold on ;
% set( gcf, 'Color', 'w', 'Units', 'Normalized', ...
%    'Position', [0.1,0.1,0.6,0.6] ) ;
% 
%    % - Plot grid.
%  plot( [xg;xg], repmat( [minlat;maxlat], 1, numel( xg )), 'Color', 0.8*[1,1,1] ) ;
%  plot( repmat( [minlon;maxlon], 1, numel( yg )), [yg;yg], 'Color', 0.8*[1,1,1] ) ;
%  
%   % - Build set of unique IDs for cells.
%  x  = minlon + 2*maxlon*rand( nx, 1 ) ;  
%  y  = minlat + 2*maxlat*rand( nx, 1 ) ;
%   xId = sum( bsxfun( @ge, x, xg(1:end-1) ), 2 ) ;
%  yId = sum( bsxfun( @ge, y, yg(1:end-1) ), 2 ) ;
%  cellId = ny * (xId - 1) + yId ;
%  
%  % - Plot cell IDs.
%  labels = arrayfun( @(k)sprintf( '%d', k ), 1:nCells, 'UniformOutput', false ) ;
%  [X,Y]  = meshgrid( (xg(1:end-1)+xg(2:end))/2, (yg(1:end-1)+yg(2:end))/2 ) ;
%  text( X(:), Y(:), labels, 'Color', 'b', 'FontSize', 8 ) ;
%  
%  
%  
%  [A B C]=textread('E.txt','%f %f %f')
%  plot( A, B, 'rx', 'LineWidth', 2, 'MarkerSize', 3 ) ;
%  
%  labels = arrayfun( @(k)sprintf( 'P%d\\in%d | %d,%d', k, cellId(k), ...
%                     C(k) ), 1:nx, 'UniformOutput', false ) ;
%                 
%                 
%   blockCount = accumarray( cellId, ones( size( cellId )), [nCells, 1] ) ;
%   blockMean_v1inRange = accumarray( cellId, C, [nCells, 1], @mean ) ;
 
% %%%%%%%%%%%%%%%%%%%
% %% Zmap format
% 
% A = xlsread('22800','A');
% A(isnan(A))=0;
% A(find(A(:,6)==0),:)=[];
% A(find(A(:,1)< 30),:)=[];
% A(find(A(:,1)> 80),:)=[];
% A(find(A(:,2)< 10),:)=[];
% A(find(A(:,2)> 45),:)=[];
% 
% A(find(A(:,6) <= 6),:)=[];
% [r c]= size(A);
% out=[];
% for i=1:r;
% %     if A(i,4) <= 1900
%     if A(i,4) >= -1250 & A(i,4) < 1800
% %      if A(i,4) >= 1964
% %      if A(i,4) >= -1250
% 
% 
% %  a =[A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),A(i,6),A(i,7),A(i,8),A(i,9)]; % Zmap
%  a =[A(i,1),A(i,2),A(i,7)];   % dpt
%  
% out = [out;a];
%     end
% end
% 
% 
% fileID = fopen('Dpt.txt','w');
% fprintf(fileID,'%f %f %f\n',out');

% fileID = fopen('22800.txt','w');
% fprintf(fileID,'%f %f %f %f %f %f %f %f %f\n',out');
% fclose(fileID);
% %%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%
% 
% 
% plot(A(600:end,4),A(600:end,7),'o')
% 
% % annual rate
% annr =[];
% for i=-1250:2006
% nn = length(find(A(:,4) == i));
% ann = [i,nn]; 
% annr = [annr;ann];
% end
% 
% plot(annr(:,1),annr(:,2))
% % average annual mag
% 
% me=[];
% for i=-1250:2006
% % [B,idx] =A(A >=i & A < i+5);
% % dd = [i+3, sum]
% [B idx]=find(A(:,4) == i);
% TF = isempty(B);
% if TF == 0;
% sb = sum(A(B(1):B(end),7))/length(B);
% mm = [i sb];
% me = [me ; mm];
% end
% end
% 
% plot(me(:,1),me(:,2))
% 
% % average annual E
% me=[];
% for i=-1250:2006
% % [B,idx] =A(A >=i & A < i+5);
% % dd = [i+3, sum]
% [B idx]=find(A(:,4) == i);
% TF = isempty(B);
% if TF == 0;
% sb = sum(1.5*(A(B(1):B(end),7))+11.8)/length(B);
% mm = [i sb];
% me = [me ; mm];
% end
% end

%%%%%%%%%%%%%%%%%%%%%%%%
% %% Zmap format/Iran
% 
% A = xlsread('IRSC','Sheet1');
% 
% [r c]= size(A);
% outA=[];
% for i=1:r;
% if 44 <= A(i,1) && A(i,1)<=62 && 24 <= A(i,2) && A(i,2) <= 40 
% a =[A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),0.922*A(i,6)+0.494,A(i,7),A(i,8),A(i,9)]; % Zmap
%  
% outA = [outA;a];
% end
% end
% 
% 
% B = xlsread ('ISC2.xlsx', 'Sheet1');
% 
% [r c]= size(B);
% outB=[];
% for i=1:r;
% 
% if 44 <= B(i,1) && B(i,1)<=62 && 24 <= B(i,2) && B(i,2) <= 40 
% if B(i,10)== 1 %mb
% 
% b =[B(i,1),B(i,2),B(i,3),B(i,4),B(i,5),1.298*B(i,6)-1.349,B(i,7),B(i,8),B(i,9)]; % Zmap
% outB = [outB;b];
% 
% elseif B(i,10)== 2 %Mw
%     
% b =[B(i,1),B(i,2),B(i,3),B(i,4),B(i,5),B(i,6),B(i,7),B(i,8),B(i,9)]; % Zmap
% outB = [outB;b];
% 
% elseif B(i,10)== 3 % ML
%     
% b =[B(i,1),B(i,2),B(i,3),B(i,4),B(i,5),0.81*B(i,6)+1.098,B(i,7),B(i,8),B(i,9)]; % Zmap
% outB = [outB;b];
% 
% elseif B(i,10)== 4 %Ms
%    
% b =[B(i,1),B(i,2),B(i,3),B(i,4),B(i,5),0.692*B(i,6)+1.945,B(i,7),B(i,8),B(i,9)]; % Zmap
% outB = [outB;b];
% 
% end
% end
% end
% 
% 
% 
% C = xlsread ('USGS.xlsx', 'Sheet1');
% 
% [r c]= size(C);
% outC=[];
% for i=1:r;
% if 44 <= C(i,1) && C(i,1)<=62 && 24 <= C(i,2) && C(i,2) <= 40 
% if C(i,10)== 1 %mb
%   
% c =[C(i,1),C(i,2),C(i,3),C(i,4),C(i,5),1.298*C(i,6)-1.349,C(i,7),C(i,8),C(i,9)]; % Zmap
% outC = [outC;c];
% 
% elseif C(i,10)== 2 %Mw
%     
% c =[C(i,1),C(i,2),C(i,3),C(i,4),C(i,5),C(i,6),C(i,7),C(i,8),C(i,9)]; % Zmap
% outC = [outC;c];
% 
% elseif C(i,10)== 3 % ML
%     
% c =[C(i,1),C(i,2),C(i,3),C(i,4),C(i,5),0.81*C(i,6)+1.098,C(i,7),C(i,8),C(i,9)]; % Zmap
% outC= [outC;c];
% 
% elseif C(i,10)== 4 %Ms
%    
% c =[C(i,1),C(i,2),C(i,3),C(i,4),C(i,5),0.692*C(i,6)+1.945,C(i,7),C(i,8),C(i,9)]; % Zmap
% outC = [outC;c];
% 
% elseif B(i,10)== 5 %lg
%    
% c =[C(i,1),C(i,2),C(i,3),C(i,4),C(i,5),0.922*C(i,6)+0.494,C(i,7),C(i,8),C(i,9)]; % Zmap
% outC = [outC;c];
% end
% end
% end
% 
% 
% [rA c]= size(outA);
% [rB c]= size(outB);
% matchId =[];
% for i=1:rA;
%  for j=1:rB;
%  
% if abs(outA(i,1)- outB(j,1)) < 0.2 && abs(outA(i,2)- outB(j,2)) < 0.2 && abs(outA(i,3)- outB(j,3))==0 && abs(outA(i,4)- outB(j,4))== 0 && (outA(i,5)-outB(j,5))==0
% %    outA(j,6)=outB(i,6)
%    matchId = [matchId;i];
% end
%  end
% end
% 
% uA = unique(matchId)
% outA(uA,:)=[];
% 
% [rC c]= size(outC);
% [rB c]= size(outB);
% matchId =[];
% for i=1:rC;
%  for j=1:rB;
%      
% if abs(outC(i,1)- outB(j,1)) < 0.2 & abs(outC(i,2)- outB(j,2)) < 0.2 & outC(i,3)== outB(j,3) & outC(i,4)== outB(j,4) & outC(i,5)== outB(j,5)
% %    outA(j,6)=outC(i,6);
%    matchId = [matchId;i],
% end
%  end
% end
% uC = unique(matchId)
% outC(uC,:)=[];
% 
% dlmwrite('orgIRSC.txt',A)
% dlmwrite('orgISC.txt',B)
% dlmwrite('orgUSGS.txt',C)
% dlmwrite('uniIRSC.txt',outA)
% dlmwrite('uniISC.txt',outB)
% dlmwrite('uniUSGS.txt',outC)
% 
% out=[outA;outB;outC];
% % fileID = fopen('IranNEW.txt','w');
% % fprintf(fileID,'%f %f %f %f %f %f %f %f %f\n','out');
% % fclose(fileID);
% 
% dlmwrite('iranUNIFIED.txt',out)
% 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% Zmap format new IRAN
% [lon lat yr mon dy mag dpt hr min typ] = textread ('ISC3.txt', '%f %f %f %f %f %f %f %f %f %s');
% A = [lon lat yr mon dy mag dpt hr min];
% fid = fopen( 'tmp.txt', 'wt');
% for ii=1:33635
%     fprintf( fid,'%f %f %f %f %f %f %f %f\n', A(ii,1), A(ii,2), A(ii,3),A(ii,4), A(ii,5), A(ii,6), A(ii,7), A(ii,8));
% end
% fclose( fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% DPT
% 
% A = xlsread('22800','A');
% A(isnan(A))=0;
% A(find(A(:,6)==0),:)=[];
% A(find(A(:,1)< 30),:)=[];
% A(find(A(:,1)> 80),:)=[];
% A(find(A(:,2)< 10),:)=[];
% A(find(A(:,2)> 45),:)=[];
% 
% [r c]= size(A);
% out=[];
% for i=1:r;
% %     if A(i,4) <= 1900
% %     if A(i,4) >= -1250 & A(i,4) < 1800
% %      if A(i,4) >= 1964
% %      if A(i,4) >= -1250
% 
%  a =[A(i,1),A(i,2),A(i,7), A(i,6)];   % dpt
%  
% out = [out;a];
% %     end
% end
% 
% 
% fileID = fopen('Dpt.txt','w');
% fprintf(fileID,'%f %f %f %f\n',out');
% fclose(fileID);

 
% fileID = fopen('iran.txt','w');
% fprintf(fileID,'%f %f %f %f %f %f %f %f %f\n',out');
% fclose(fileID);
 




% %%%%%%%%%%%%%%%%%%%%%%
% [lat lon exx eyy exy vorcit RL LL e1 e2 azi] = textread ('E.txt', '%f %f %f %f %f %f %f %f %f %f %f');
% A = [lon lat e1 e2];
% r = length(A);
% 
% out=[];
% for i=1:r;
%     bb=i
%     inv = sqrt(A(i,3)^2 + A(i,4)^2);
%     st = (A(i,3)+A(i,4))/max(abs(A(i,3)), abs(A(i,4)));
%     pr =[A(i,1) A(i,2) inv st];
%     out=[out;pr];
% end 
% 
% fileID = fopen('EE.txt','w');
% fprintf(fileID,'%f %f %f %f\n',out'); %D
% fclose(fileID);
% 
% 
% 
% % %%%%%%%%%%%%
% [lat lon b] = textread ('160_b.txt', '%f %f %f');
% A = [lon lat b];
% % b(isnan(b))=-100;%
% A(find(A(:,3) == -100),:)=[];
% A(isnan(A(:,3)),:)=[];
%  
% fileID = fopen('bval.txt','w');
% fprintf(fileID,'%f %f %f\n',A'); %D
% fclose(fileID);

% [lon lat yr mon dy mag dpt hr min] = textread ('IRAN(unified+KUdeclust).txt', '%f %f %f %f %f %f %f %f %f');
% m = 1.5*(mag)+9.1;
% out = [lon lat m];
% 
% o=[];
% % %  for d=30:0.5:80
% % %     for mm = 10:0.5:45
% % %      oo = [d, mm, 0];
% % %      out = [out; oo];
% % %     end
% % %  end
% % 
% % oo=[];
% % 
% % for i=1:r;
% % 
% % % if ((A(i,3) >= -1250) && (A(i,3) < 1900) && (A(i,6) >= 5.0))
% % %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% % %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % % out = [out;a];
% % % oo = [oo;b];
% %  
% % if ((A(i,3) >= 1900) && (A(i,3) < 1927) && (A(i,6) >= 5.0))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% % elseif ((A(i,3) >= 1927) && (A(i,3) < 1963) && (A(i,6) >= 4.5))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% %  
% % elseif ((A(i,3) >= 1963) && (A(i,3) < 1973) && (A(i,6) >= 4.8))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% % elseif ((A(i,3) >= 1973) && (A(i,3) < 1995) && (A(i,6) >= 4.5))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% % elseif  ((A(i,3) >= 1995) && (A(i,3) < 2001) && (A(i,6) >= 4.2))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% % elseif ((A(i,3) >= 2001) && (A(i,3) < 2007) && (A(i,6) >= 4.0))
% %  a =[A(i,1),A(i,2),(1.5*A(i,6)+11.8)]; % E
% %  b =[A(i,1),A(i,2),A(i,3)];   %year
% % out = [out;a];
% % oo = [oo;b];
% % 
% %     end
% % end
% % 
% % 
% % % %% magnitude /time distribution
% % % plot(oo(:,3),oo(:,1),'o')
% % % grid on
% % % xlabel('Time (years)','FontSize',18,'FontWeight','bold')
% % % ylabel('Longitude (degree)','FontSize',18,'FontWeight','bold')
% % % set(gca,'FontSize',15)
% % % 
% % % % xlim([1995,2001])
% % % % title('4.8 - 5.4')
% % 
% % %   mm =sum(out(:,3))/length(out);
% % %  out(:,3) = out(:,3)-mm;
% % %  mino=min(out(:,3));
% % %  
% % %  for d=30:0.5:80
% % %     for mm = 10:0.5:45
% % %      oo = [d, mm, mino];
% % %      out = [out; oo];
% % %     end
% % %  end
% % 
% %  
% % delete E.txt
% % fileID = fopen('E.txt','w');
% % % fprintf(fileID,'%f %f %f %f\n',out'); %E
% % fprintf(fileID,'%f %f %f \n',out'); %D
% % fclose(fileID);
% % % movefile('E.txt','../eng3/.') 
% % % system('./eng3/runEng.sh')
% 
% fileID = fopen('E.txt','w');
% fprintf(fileID,'%f %f %f\n',out');
