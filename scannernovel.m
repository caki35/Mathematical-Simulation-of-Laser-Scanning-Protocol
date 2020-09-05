function [N,maksvalue,minvalue]= scannernovel(wait_duration,fscan,freprate,small_location_derivation,higher_location_derivation,timedifpoint,stdwm)
%----- sine signal is added on to a ramp signal -----%
%-------------------------- input arguments ------------------------%
%       wait_duration: how long scanning will continue              %
%       fscan: frequency of the sine wave                           %
%       freprate: frequency of the ramp wave                        %
%       small_location_derivation = amplitude of the sine signal    %
%       higher_location_derivation = amplitude of the ramp signal   %
%       timedifpoint: the point at which time dif are calculated    %
%       stdwm = show time differences without modulation(if 1)      %
%-------------------------------------------------------------------%
%-------------------------------------------------------------------%

fs = 40000*fscan;                                           % sampling frequency
dt=1/fs;                                                    % time interval
t = 0:dt:wait_duration;                                     % time vector

%the function for scannig
%it will keep amplitude voltage of scanner between 0V-ramp signal
A1=small_location_derivation;
A2=higher_location_derivation;
scanner=A1*sin(2*pi*t*fscan);
k=(A2-2*A1)/A2;
rampscan=((A2+A2*sawtooth(2*pi*t*freprate))/2)*k+A1;
%the signal which we try to obtain ultimately is refered as y through the code 
y=rampscan+scanner; 
rampscan_wo_modulation=A2*((sawtooth(2*pi*t*freprate)+1)/2);


%% Calculate the change in time difference between consecutive points 
% Returns Zero-Crossing Indices Of Argument Vector
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); 
% Returns Indices Of Argument Vector at the specifict point
% timedifpoint is taken from user
zt = zci(y-timedifpoint);  


% the times at which the scanner cross the timedifpoint
t_crosszero=t(zt');          
%While t_crosszero will be used up-down scan backup is used to up and down
t_crosszero_backup=t_crosszero;

delta_tz=zeros(length(t_crosszero)-1,1);          % initiate time differences   
for i=1:length(delta_tz)                          % time differences calculation
delta_tz(i)=t_crosszero(i+1)-t_crosszero(i);
end
N=length(t_crosszero(1:2:end));                             % data points

%Note length(delta_tz)=length(t_crosszero)-1
% remove firt cross. There is no time difference
t_crosszero(1)=[];

%Those variable will be used to draw time dif while consider both up-down scan
t_dif1=t_crosszero;
d_dif1=delta_tz;

%detect outliers which are more bigger
%our liers consist of time difference between.. 
%..the last point of a scan and the first point of another scan
u=mean(d_dif1);
var=std(d_dif1);
index_to_remove=find(abs(d_dif1-u)>3*var);
t_dif1(index_to_remove)=[];
d_dif1(index_to_remove)=[];

maksvalue=max(d_dif1);
minvalue=min(d_dif1);

%Those variable will be used to draw time dif while consider up scan
t_dif2=t_crosszero_backup(1:2:end);

d_dif2=zeros(length(t_dif2)-1,1);          % initiate time differences   
for i=1:length(d_dif2)                     % time differences calculation
d_dif2(i)=t_dif2(i+1)-t_dif2(i);
end

% remove firt cross. There is no time difference
t_dif2(1)=[];

%detect outliers which are more bigger
u=mean(d_dif2);
var=std(d_dif2);
index_to_remove=find(abs(d_dif2-u)>3*var);
t_dif2(index_to_remove)=[];
d_dif2(index_to_remove)=[];

%Those variable will be used to draw time dif while consider down scan
t_dif3=t_crosszero_backup(2:2:end);

d_dif3=zeros(length(t_dif3)-1,1);          % initiate time differences   
for i=1:length(d_dif3)                     % time differences calculation
d_dif3(i)=t_dif3(i+1)-t_dif3(i);
end

% remove firt cross. There is no time difference
t_dif3(1)=[];


u=mean(d_dif3);
var=std(d_dif3);
index_to_remove=find(abs(d_dif3-u)>3*var);
t_dif3(index_to_remove)=[];
d_dif3(index_to_remove)=[];

%do you want show time dif
if stdwm==1
    % timedifpoint is taken from user
    zt2 = zci(rampscan_wo_modulation-timedifpoint);  %here ramp_without_sine is used instead of y
    % the times at which the scanner cross the timedifpoint
    t_crosszero_withoutsine=t(zt2');

    delta_tz2=zeros(length(t_crosszero_withoutsine)-1,1);          % initiate time differences   
    for i=1:length(delta_tz2)                          % time differences calculation
    delta_tz2(i)=t_crosszero_withoutsine(i+1)-t_crosszero_withoutsine(i);
    end
    % remove firt cross. There is no time difference
    t_crosszero_withoutsine(1)=[];
    t_crosszero_withoutsine(1:2:end)=[];        % remove negative gradient zero cross(return zero)
    delta_tz2(1:2:end)=[];           % remove time differences at which negative gradient zero cross occured(return zero)
end;

%% Plotting first figure
figure(1)
subplot(3,2,1)
plot(t,scanner)
title('Sine Signal')
xlabel('Time (sec.)'),ylabel('Position(Sine)');
grid on

subplot(3,2,3)
plot(t,rampscan)
title('Ramp Signal')
xlabel('Time (sec.)'),ylabel('Position(Ramp)');
grid on

subplot(3,2,5)
plot(t,y);
title('Ramp+Sine Signal')
xlabel('Time (sec.)'), ylabel('Position');
grid on

%Plot time differences
subplot(3,2,2)
scatter(t_dif1,d_dif1,'MarkerEdgeColor',[0 0 0.7],'MarkerFaceColor',[0.1 0.1 1])
title('Time Differences for Bidirectional Scanning')
xlabel('Time (sec.)'), ylabel('Time Dif (sec.)');
grid on

subplot(3,2,4)
scatter(t_dif2,d_dif2,'MarkerEdgeColor',[0.7 0 0],'MarkerFaceColor',[1 0 0])
title('Time Differences for Unidirectional Scanning (Up Scanning)')
xlabel('Time (sec.)'), ylabel('Time (sec.)');
grid on

subplot(3,2,6)
scatter(t_dif3,d_dif3,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7])
title('Time Differences for Unidirectional Scanning (Down Scanning)')
xlabel('Time (sec.)'), ylabel('Time (sec.)');
grid on

%The sine wave without ramp
figure(2)
subplot(1,5,1)
plot(t,scanner)
xlabel('Time (sec.)'),ylabel('Position(Sine)');
camroll(270)

%The ramp signal without sine
subplot(1,5,2)
plot(t,rampscan)
xlabel('Time (sec.)'),ylabel('Position(Ramp)');
camroll(270)

%Rampe+sine
subplot(1,5,3)
plot(t,y)
xlabel('Time (sec.)'),ylabel('Position(Ramp+Sine)');
grid on
camroll(270)



%% The Last adding
%plot time differences for:
    %going
    %returning
    %both
%Adding ability to select the point at which time differences is calculated
figure(3)

%plot ramp+sine signal
subplot(1,3,1)
plot(t,y)
xlabel('Time (sec.)'),ylabel('Position(Ramp+Sine)');
grid on
camroll(270)

%Plot time-dif in up+down scan
subplot(1,3,2)
scatter(t_dif1,d_dif1)
xlabel('Time (sec.)'), ylabel('Time Differences(sec.) for B.S.');
grid on
camroll(270)

%Plot time-dif in up and down scan on same plot with different color
subplot(1,3,3)
scatter(t_dif2,d_dif2)
hold on 
scatter(t_dif3,d_dif3)
if stdwm==1
    hold on
    scatter(t_crosszero_withoutsine,delta_tz2,'MarkerEdgeColor','g')
end;
xlabel('Time (sec.)'), ylabel('Time Differences(sec.) for U.S.');
grid on
camroll(270)


