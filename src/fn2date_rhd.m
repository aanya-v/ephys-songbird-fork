function [hour,day,month,year]=fn2date(fn);
%[hour,day,month,year]=fn2date(fn);
%

p = findstr(fn,'.rhd');
if (length(p)<1)
    p=findstr(fn,'.mat');
end
if (length(p)<1)
    p=findstr(fn,'.wav');
end

% if (length(p)<1)
%     disp(['not rhd, mat or wav file?']);
%     return;
% end
% 
% p2=findstr(fn,'.');
% p3=find(p2==p(end)); 
% if (length(p3)<1)|(length(p2)<2) %if fn doesn't include '.rhd'
%     disp(['weird fn = ',fn]);
%     return;
% end
% 
% p = p2(p3-1);


% t = fn([(p(end)-6):(p(end)-1)]);
% dt = fn([(p(end)-13):(p(end)-8)]);

idx = strfind(fn,'_'); 
dt = fn(idx(1)+1:idx(2)-1); %after first "_" is the date
t = fn(idx(2)+1:idx(3)-1);

hour = str2num(t(1:2)) + str2num(t(3:4))/60.0;
day   = str2num(dt(5:6));
month = str2num(dt(3:4));
year  = 2000 + str2num(dt(1:2));
return;
