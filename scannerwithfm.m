function [N,maksvalue,minvalue]=scannerwithfm(wait_duration,fc,fm,delta_f,timedifpoint,stdwm)
%----------------------- frequency modulation ----------------------%
%-------------------------- input arguments ------------------------%
%       wait_duration: how long scanning will continue              %
%       fc: central frequency of the modulated signal               %
%       fm: frequency of modulation                                 %
%       delta_f : modulation bandwidth                              %
%       timedifpoint: the point at which time dif are calculated    %
%       stdwm = show time differences without modulation(if 1)      %
%-------------------------------------------------------------------%
%-------------------------------------------------------------------%

fs = 50000*fc;                                                              % sampling frequency
t = 0:1/fs:wait_duration;                                                   % time vector
ramp_original=sawtooth(2*pi*fc*t,1/2);                                          % original signal which will be modulated
modulation_index=delta_f/fm;                                                % modulation index
ramp_modulated = sawtooth((2*pi*fc*t)+(modulation_index.*sin(2*pi*fm*t)),1/2);  % modulation 
[pks,locs] = findpeaks(ramp_modulated);                                     % find peaks
N=length(pks);                                                              % how many sweeps occured


%% Calculate the change in time difference between consecutive points 
% Returns Zero-Crossing Indices Of Argument Vector
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); 
% Returns Indices Of Argument Vector at the specifict point
% timedifpoint is taken from user
zt = zci(ramp_modulated-timedifpoint);  

% the times at which the scanner cross the timedifpoint
t_crosszero=t(zt');     

%While t_crosszero will be used up-down scan backup is used to up and down
t_crosszero_backup=t_crosszero;

delta_tz=zeros(length(t_crosszero)-1,1);          % initiate time differences   
for i=1:length(delta_tz)                          % time differences calculation
delta_tz(i)=t_crosszero(i+1)-t_crosszero(i);
end

%Note length(delta_tz)=length(t_crosszero)-1
% remove firt cross. There is no time difference
t_crosszero(1)=[];

%Those variable will be used to draw time dif while consider both up-down scan
t_dif1=t_crosszero;
d_dif1=delta_tz;

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


%Those variable will be used to draw time dif while consider down scan
t_dif3=t_crosszero_backup(2:2:end);

d_dif3=zeros(length(t_dif3)-1,1);          % initiate time differences   
for i=1:length(d_dif3)                     % time differences calculation
d_dif3(i)=t_dif3(i+1)-t_dif3(i);
end

% remove firt cross. There is no time difference
t_dif3(1)=[];

%do you want show time dif
if stdwm==1
    % timedifpoint is taken from user
    zt2 = zci(ramp_original-timedifpoint);  %here ramp_without_sine is used instead of y
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

%% Frequency Analysis of Modulated signal
% vector of frequencies in Hz
npnts = wait_duration*fs+1;
hz = linspace(0,fs/2,floor(npnts/2)+1);

% amplitude spectrum via Fourier transform
Signal = ramp_modulated-mean(ramp_modulated);
power_spectrum = abs(fft(Signal)/npnts).^2;

% now for a time-frequency analysis
% use MATLAB's spectrogram function (in the signal processing toolbox)
pks_delta_t = findpeaks(d_dif1);                                           % find peaks 
sizeratio=length(pks_delta_t)*5;
nsc = floor(length(Signal)/sizeratio);                                      % Hanning window size
nov = floor(nsc/2);                                                         % Overlapping between adajacent windows
nff = max(256,2^nextpow2(nsc));                                             % number of samples to compute the FFT
[powspect,frex,time] = spectrogram(Signal,nsc,nov,nff,fs);


figure(1)
subplot(3,2,1)
plot(t,ramp_original);
title('original signal without modulation')
xlabel('Time (sec.)'), ylabel('Amplitude');
grid

subplot(3,2,3)
plot(t,sin(2*pi*fm*t));
title('modulation signal')
xlabel('Time (sec.)'), ylabel('Amplitude');
grid
 
subplot(3,2,5)
plot(t,ramp_modulated);
title('modulated signal')
xlabel('Time (sec.)'), ylabel('Amplitude');
grid 

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

figure(2)
subplot (1,2,1); 
plot(hz,power_spectrum(1:length(hz)),'linew',2)
set(gca,'xlim',[0 (fc+delta_f)*5])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Power Spectrum')


% show the time-frequency power plot
subplot(1,2,2)
imagesc(time,frex,abs(powspect).^2)
axis xy
set(gca,'ylim',frex([1 dsearchn(frex,(fc+delta_f)*5)]),'xlim',time([1 end]));
xlabel('Time (sec.)'), ylabel('Frequency (Hz)')
title('Spectrogram');
colormap hot



%% The Last adding
%plot time differences for:
    %going
    %returning
    %both
%Adding ability to select the point at which time differences is calculated
figure(3)

%plot ramp+sine signal
subplot(1,3,1)
plot(t,ramp_modulated);
xlabel('Time (sec.)'), ylabel('Position');
grid 
camroll(270)

%Plot time-dif in up+down scan
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
if stdwm==1
    hold on
    scatter(t_crosszero_original,delta_tz2,'MarkerEdgeColor','g')
end;
xlabel('Time (sec.)'), ylabel('Time Differences(sec.) for U.S.');
legend('Time Differences In Up-Scan','Time Differences In Down-Scan')
grid on
camroll(270)
