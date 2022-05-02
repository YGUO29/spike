%% you need to change most of the paths in this block

addpath(genpath('D:\MatlabToolbox\Kilosort-main_3.0')) % path to kilosort folder
addpath('D:\MatlabToolbox\npy-matlab-master') % for converting to Phy
rootZ = 'D:\yueqi\axoft-data\cortexlab-dataset'; % the raw data binary file is in this folder
rootH = 'D:\yueqi\axoft-data\cortexlab-dataset'; % path to temporary binary file (same size as data, should be on fast SSD)
pathToYourConfigFile = 'D:\MatlabToolbox\Kilosort-main_3.0\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = 'neuropixPhase3A_kilosortChanMap.mat';

ops.trange    = [0 Inf]; % time range to sort
ops.NchanTOT  = 385; % total number of channels in your recording

run(fullfile(pathToYourConfigFile, 'configFile384.m'))
ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);

%% load raw waveform
fn = fullfile(rootH, 'rawDataSample.bin');
fid = fopen(fn, 'r');
dat = fread(fid, [385 Inf], '*int16');
fclose(fid);

chanMap = readNPY(fullfile(rootH, 'channel_map.npy'));
dat = dat(chanMap+1,:);

% visualize raw waveform
figure('color','w'); 
fs = 3e4;
ch_idx = 1:10:374;
t = [1:30000]./fs;
for i = 1:length(ch_idx)
    hold on, plot(t, dat(ch_idx(i),1:30000)' + 60*(i-1))
end
xlabel('time (s)')
ylabel('channel number')
set(gca, 'ytick', 60.*[1:length(ch_idx)],...
    'yticklabels', arrayfun(@num2str, ch_idx, 'UniformOutput', false))
ylim([-60, 2300])
% imagesc(dat(:,1:30000)), colormap(gray)

%% this block runs all the steps of Kilosort
fprintf('Looking for data inside %s \n', rootZ)

% main parameter changes from Kilosort2 to v2.5
ops.sig        = 20;  % spatial smoothness constant for registration
ops.fshigh     = 300; % high-pass more aggresively
ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 

% main parameter changes from Kilosort2.5 to v3.0
ops.Th       = [9 9];

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
% ops.fbinary = fullfile(rootZ, fs(1).name);
ops.fbinary        = fs(1).name;
ops.root           = rootZ;

rez                = preprocessDataSub(ops);
rez                = datashift2(rez, 1);

[rez, st3, tF]     = extract_spikes(rez);

rez                = template_learning(rez, tF, st3);

[rez, st3, tF]     = trackAndSort(rez);

rez                = final_clustering(rez, tF, st3);

rez                = find_merges(rez, 1);

rootZ = fullfile(rootZ, 'kilosort3');
mkdir(rootZ)
rezToPhy2(rez, rootZ);

%% load the results after clustering (.py data) into matlab
addpath(genpath('D:\MatlabToolbox\spikes-master'))
myKsDir = 'D:\yueqi\axoft-data\cortexlab-dataset\kilosort3';
sp = loadKSdir(myKsDir);

%% select a cluster, plot raw waveform
iCluster = 330;
spike_time = sp.st(sp.clu == iCluster-1);
spike_window = [-40 30];

dat_chunk = zeros(sp.n_channels_dat, spike_window(2)-spike_window(1)+1, length(spike_time));
for iSpike = 1:length(spike_time)
    dat_chunk(:,:, iSpike) = double(dat(:, ...
        round(sp.sample_rate*spike_time(iSpike))+spike_window(1) : ...
        round(sp.sample_rate*spike_time(iSpike))+spike_window(2)));
end
dat_mean = squeeze(mean(dat_chunk,3));
[~, ind_ch] = max(peak2peak(dat_mean,2));

% visualize the mean waveform on all channels
t = [1:size(dat_mean,2)]./fs.*1000;
figure('color','w')
imagesc(t, 1:size(dat_mean,1), dat_mean, [-max(abs(dat_mean(:))), max(abs(dat_mean(:)))]);
colormap(jet), colorbar
set(gca,'Ydir','normal')
title(['cluster = ', num2str(iCluster), ', quality = ', num2str(sp.cgs(iCluster))])
xlabel('time (ms)')
ylabel('channel number')

% visualize the waveform
figure('color','w')
for iSpike = 1:length(spike_time)
    hold on, plot(t, dat_chunk(ind_ch, :, iSpike))
end
hold on, plot(t, dat_mean(ind_ch, :),'linewidth',2,'color','k')
title(['cluster = ', num2str(iCluster), ', quality = ', num2str(sp.cgs(iCluster))])
xlabel('time (ms)')
ylabel('Voltage (uV)')
%% raster plot
xx = []; yy = [];
for iCluster = 1:max(sp.clu+1)
    iCluster = double(iCluster);
    yy = [yy; iCluster.*ones(length(find(sp.clu == iCluster-1)), 1)];
    xx = [xx; sp.st(sp.clu == iCluster-1)];
end

figure('color','w')
scatter(xx, yy, '.')
xlabel('time (s)')
ylabel('Cluster number')
xlim([0 30]), ylim([1 max(yy)])

%% psth
clu_idx = 400:500;
sp_train = sp.st(sp.clu >=clu_idx(1) & sp.clu <=clu_idx(end));
% sp_train = sp.st;
nCluster = length(clu_idx);
window = 1e-2; % bin width in seconds: 10ms
tend = [];

psth = spike2psth(sp_train, window, nCluster, tend);
xlabel('time (s)')
ylabel('number of spikes/second')
xlim([0 30])
title(['cluster ', num2str(clu_idx(1)), '-', num2str(clu_idx(end))])
%% auto- and crosscorrelogram
clu_number = [2 3];

delay = -0.05:1e-3:0.05; % seconds
sp_train1 = sp.st(sp.clu == clu_number(1));
sp_train2 = sp.st(sp.clu == clu_number(2));
window = 1e-3; % 1ms window
tend = [];
psth1 = spike2psth(sp_train1, window, 1, tend);
psth2 = spike2psth(sp_train2, window, 1, tend);
[r, lags] = xcorr(psth1, psth2, 50);
r = r./sum(r);
if clu_number(2) == clu_number(1)
    r(51) = 0;
end

figure('color','w')
subplot(2,1,1)
scatter(sp_train1, ones(1,length(sp_train1)),'.')
hold on
scatter(sp_train2, 2.*ones(1,length(sp_train2)),'.')
ylim([0 3])
xlabel('time (s)')
ylabel('2 clusters')
title('raster plot')
subplot(2,1,2)
bar(lags, r, 'FaceColor', 'k')
xlabel('lag (ms)')
if clu_number(2) == clu_number(1)
    title(['auto-correlogram, cluster number = ', num2str(clu_number(1)), ...
        ', quality = ', num2str(sp.cgs(clu_number(1)+1))])
else
    title(['cross-correlogram, cluster number = ', num2str(clu_number),...
        ', quality = ', num2str(sp.cgs(clu_number+1))])
end


