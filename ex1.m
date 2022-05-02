% exercise 1: nsx files

addpath(genpath('D:\MatlabToolbox\NPMK-5.5.2.0'))
addpath(genpath(cd))

para.datapath = 'D:\yueqi\axoft-data';
para.filename = {'12_072521_002.ns6', '12_072521_002.ns4'}; 
para.fs = [3e4, 1e4];
openNSx('report', 'read',fullfile(para.datapath, para.filename{1}));
openNSx('report', 'read',fullfile(para.datapath, para.filename{2}));

%% plot waveform
% raw data
t = 1/para.fs(1):1/para.fs(1):size(NS6.Data,2)/para.fs(1);

figure('color', 'w')
subplot(2,2,1)
plot(t, double(NS6.Data(1,:)).*NS6.ElectrodesInfo(1).Resolution)
ylim([-350 350]), xlim([0 t(end)])
title('channel 1, raw data')
xlabel('time (s)')
ylabel('voltage (uV)')
subplot(2,2,2)
plot(t, double(NS6.Data(32,:)).*NS6.ElectrodesInfo(1).Resolution)
ylim([-350 350]), xlim([0 t(end)])
title('channel 32, raw data')
xlabel('time (s)')
ylabel('voltage (uV)')

% filtered data
t = 1/para.fs(2):1/para.fs(2):size(NS4.Data,2)/para.fs(2);

subplot(2,2,3)
plot(t, double(NS4.Data(1,:)).*NS4.ElectrodesInfo(1).Resolution)
ylim([-200 200]), xlim([0 t(end)])
title('channel 1, filtered data')
xlabel('time (s)')
ylabel('voltage (uV)')

subplot(2,2,4)
plot(t, double(NS4.Data(32,:)).*NS4.ElectrodesInfo(1).Resolution)
ylim([-200 200]), xlim([0 t(end)])
title('channel 32, filtered data')
xlabel('time (s)')
ylabel('voltage (uV)')


%% spike detection on NS4

t = 1/para.fs(2):1/para.fs(2):size(NS4.Data,2)/para.fs(2);
raw_wav = double(NS4.Data(32,:));

% set a threshold (below baseline to detect negative voltage changes)
threshold = mean(raw_wav) - 8*std(raw_wav);
timestamps = edgeDetect(raw_wav, threshold, 'falling'); % spike time 
timestamps_s = timestamps./para.fs(2);

figure('color','w'),
subplot(2,1,1)
plot(t, raw_wav), hold on
scatter(timestamps_s, threshold.*ones(1,length(timestamps)),'.')
xlabel('time (s)')
ylabel('voltage (uV)')
title('filtered waveform')
xlim([0 max(t)])

subplot(2,1,2)
scatter(timestamps_s, ones(1,length(timestamps)),50, '.')
xlabel('time (s)')
title('raster plot')
xlim([0 max(t)])

% psth = spike2psth(timestamps, 1e-2, 600);

% sanity check
% check spike waveforms
nSpikes = length(timestamps);
clear wav_segment
figure('color','w'),
% extract 3ms before and after each timestamp
cutoff = round(3e-3*para.fs(2));
t = [1:2*cutoff+1]./para.fs(2).*1000;
for i = 1:nSpikes
    wav_segment(i,:) = raw_wav(timestamps(i)-cutoff :...
        timestamps(i)+cutoff);
    hold on, plot(t, wav_segment(i,:))
end
hold on, plot(t, mean(wav_segment),'linewidth',2,'color','k')
xlabel('time (ms)')
ylabel('voltage (uV)')
title('raw waveforms of detected spikes')

% check inter-spike interval
figure('color','w'),
isi = diff(timestamps)./para.fs(2).*1e3;
histogram(isi(isi<10),100)
title('inter-spike-intervals')
xlabel('time (ms)')
min(isi) % if recorded single unit, this should be >2ms

