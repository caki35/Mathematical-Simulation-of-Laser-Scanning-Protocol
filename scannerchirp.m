function [N,maksvalue,minvalue]= scannerchirp(wait_duration,fsweep,f_initial,f_final,duty_cycle_percent,timedifpoint,stdwm)
%------------------------ frequency sweeping -----------------------%
%-------------------------- input arguments ------------------------%
%       wait_duration: how long scanning will continue              %
%       fsweep: frequency of original signal - sine wave            %
%       f_initial: initial frequency of sweeping                    %
%       f_final: initial frequency of sweeping                      %
%       duty_cycle_percent:                                         %
%       timedifpoint: the point at which time dif are calculated    %
%       stdwm = show time differences without modulation(if 1)      %
%-------------------------------------------------------------------%
%--------------- Report file name: Chirp report --------------------%
%-------------------------------------------------------------------%
%% Parameters
%sweep from initial frequency to finial frequency
B=f_final-f_initial;                                        % Bandwidth
f_center=(f_final+f_initial)/2;                             % center frequency
fs = 4000*f_final;                                          % sampling frequency
dt=1/fs;                                                    % time interval
t = 0:dt:wait_duration;                                     % time vector
Tsweep=1/fsweep;                                            % sweep time
duty_cycle=duty_cycle_percent/100;                          % duty cycle
T_repetition=Tsweep/duty_cycle;                             % Duty Cycle=Tsweep/Trepetition
d=0:T_repetition:wait_duration;                             % delay vector for pulse train
alfa=B/Tsweep;                                              % sweep rate


%% Creating signals
unitstep1=(t>=0);                                          %u(t)
unitstep2=(t-Tsweep>0);                                    %u(t-T)
rectangular=unitstep1-unitstep2;                           %Rect(t/T)=u(t)-u(t-T)
f_instant=rectangular.*(alfa*t+f_initial);                 %instant frequency              
%instant frequency with fsweep repetition rate
f_instant= pulstran(t,d,f_instant,fs);

% chirped signal from directly formula
chirp_signal=rectangular.*sin(2*pi*(f_initial*t+(t.^2)*B/(2*Tsweep))); 
%pulse train with chirp signal
chirp_signal_train = pulstran(t,d,chirp_signal,fs);
% signal oscilating at center frequency without chirp
original_signal=cos(2*pi*f_center*t);
original_signal_onescan=rectangular.*original_signal;
original_signal_train=pulstran(t,d,original_signal_onescan,fs);

%take only one window to obtain sweep for powerspectrum
a=nonzeros(rectangular)';
chirp_signal_window=chirp_signal(1:length(a));

%% Calculate the change in time difference between consecutive points 
% Returns Zero-Crossing Indices Of Argument Vector
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); 
% Returns Indices Of Argument Vector at the specifict point
% timedifpoint is taken from user
zt = zci(chirp_signal_train-timedifpoint);  


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
    zt2 = zci(original_signal_train-timedifpoint);  %here ramp_without_sine is used instead of y
    % the times at which the scanner cross the timedifpoint
    t_crosszero_original=t(zt2');

    delta_tz2=zeros(length(t_crosszero_original)-1,1);          % initiate time differences   
    for i=1:length(delta_tz2)                          % time differences calculation
    delta_tz2(i)=t_crosszero_original(i+1)-t_crosszero_original(i);
    end
    % remove firt cross. There is no time difference
    t_crosszero_original(1)=[];
    %detect outliers which are more bigger
    %our liers consist of time difference between.. 
    %..the last point of a scan and the first point of another scan
    u=mean(delta_tz2);
    var=std(delta_tz2);
    index_to_remove=find(abs(delta_tz2-u)>3*var);
    t_crosszero_original(index_to_remove)=[];
    delta_tz2(index_to_remove)=[];
end;


%% Power spectrum 
% vector of frequencies in Hz
npnts = t(length(a))*fs+1;                     %time points
hz = linspace(0,fs/2,floor(npnts/2)+1);        %frequency point
% amplitude spectrum via Fourier transform
Signal = detrend(chirp_signal_window);         %normalization
power_spectrum = abs(fft(Signal)/npnts).^2;    %power spectrum



% use MATLAB's spectrogram function (in the signal proacessing toolbox)
peaks_deltat= findpeaks(delta_tz);       %find peaks
sizeratio=length(peaks_deltat)*10;        %Signal-window size ratio
nsc = floor(length(chirp_signal_train)/sizeratio);  %Hanning window size
nov = floor(nsc/2);                     %Set overlapping between adjacent windows %50 
nff = max(256,2^nextpow2(nsc));         %number of samples to compute the FFT
[powspect,frex,time] = spectrogram(chirp_signal_train,nsc,nov,nff,fs);

%% Ploting
subplot(3,2,1)
plot(t,original_signal_train);
title('Original Signal')
xlabel('Time (sec.)'), ylabel('Position');
grid on

subplot(3,2,3)
plot(t,f_instant);
title('Frequency Sweeping')
xlabel('Time (sec.)'), ylabel('Frquency');
grid on

subplot(3,2,5)
plot(t,chirp_signal_train);
title('Chirp Signal')
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

%% show the power spectrum
figure(2)
subplot (2,1,1); 
plot(hz,power_spectrum(1:length(hz)),'-*','MarkerFaceColor',[1 .6 .6],'MarkerSize',3,'MarkerEdgeColor','red','linew',2);
set(gca,'xlim',[0 (f_final)*3]);clc
xlabel('Frequency (Hz)'), ylabel('Power');
title('Power Spectrum');

% %  show the spectogram
subplot (2,1,2); 
imagesc(time,frex,abs(powspect).^2);
axis xy
set(gca,'ylim',frex([1 dsearchn(frex,(f_final)*3)]),'xlim',time([1 end]));
xlabel('Time (sec.)'), ylabel('Frequency (Hz)');
title('Spectrogram');
colormap hot;



%% The Last adding
%plot time differences for:
    %going
    %returning
    %both
%Adding ability to select the point at which time differences is calculated
figure(3)

%plot ramp+sine signal
subplot(1,3,1)
plot(t,chirp_signal_train)
xlabel('Time (sec.)'),ylabel('Position(Chirp)');
grid on
camroll(270)

%Plot time-dif in double scan
subplot(1,3,2)
scatter(t_dif1,d_dif1)
xlabel('Time (sec.)'), ylabel('Time Differences(sec.) for B.S.');
grid on
camroll(270)


%Plot time-dif in up+down scan on same plot with different color
subplot(1,3,3)
scatter(t_dif2,d_dif2)
hold on 
scatter(t_dif3,d_dif3)
xlabel('Time (sec.)'), ylabel('Time Differences(sec.) for U.S.');
if stdwm==1
    hold on
    scatter(t_crosszero_original,delta_tz2,'MarkerEdgeColor','g')
end;
grid on
legend('Time Differences In Up-Scan','Time Differences In Down-Scan')
camroll(270)

