clc
clear all
close all
addpath(genpath('projects_common_material/biosig'));
filename1 = 'project2-data-example\anonymous.20170613.161402.offline.mi.mi_bhbf.gdf';
filename2 = 'project2-data-example\anonymous.20170613.162331.offline.mi.mi_bhbf.gdf';
filename3 = 'project2-data-example\anonymous.20170613.162934.offline.mi.mi_bhbf.gdf';
[s1, h1] = sload(filename1);
[s2, h2] = sload(filename2);
[s3, h3] = sload(filename3);

s = [s1;s2;s3];
size_ = [length(s1),length(s2),length(s3)];
label = [ones(size_(1),1); ones(size_(2),1)*2; ones(size_(3),1)*3];
s = [label, s];

h2.EVENT.POS = h2.EVENT.POS + length(s1);
h3.EVENT.POS = h3.EVENT.POS + length(s2) + length(s1);

%% Isolate interesting bandpass

%faire un filtre et isole les bandes mu et beta
sample_rate = h1.SampleRate;

% Bande mu = [7-14] Hz
[b,a] = butter(4,[7,14]/sample_rate/2);
fvtool(b,a)
mu_data = filter(b,a,s(:,2:end));

% Bande beta = [15-35] Hz
[b,a] = butter(4,[15,35]/sample_rate/2);
fvtool(b,a)
beta_data = filter(b,a,s(:,2:end));

%% Spatial filtering:

% CAR filtering:
car_mu_data = [];
for i = 1:length(mu_data)
    moy = mean(mu_data(i,1:end-1));
    for j = 1:size(mu_data,2)
        car_mu_data(i,j) = mu_data(i,j) - moy;
    end
end

%% Laplacian filter
 
% Small Laplacian filtering
small_laplacian = load('laplacian_16_10-20_mi.mat');
small_lap_mu_data = mu_data(:,1:end-1) * small_laplacian.lap;

% Large laplacian filter
large_laplacian = load('big_laplacian.mat');
large_lap_mu_data = mu_data(:,1:end-1) * small_laplacian.lap;

%% Power spectral density estimate

psd = pwelch(s(:,9),256,128,[],512)
figure()
plot(log(psd))
% we see alpha peak
% We must check it looks like this for our own signals.

% We get PSD for all temporal signal: but we don't want that:
% We need a temporal window for FFT
% Buffer filled every 62.5 ms --> compute spatial filter and then PSD every 62.5 ms
% This corresponds to 60 Hz. To do on 1s.
% stop position = start position + 1s.
%Fixation = 786
%Both Hands/Right Cue = 773;
%Both Feet/Left Cue = 771;
%Continuous Feedback = 781;
%Boom Target Hit = 897;
%Boom Target Miss = 898;


EventId = 773;
StartPositions = h1.EVENT.POS(h1.EVENT.TYP==EventId);
StopPositions = StartPositions + 511;
psd = [];

NumChannels = size(s,2)-2;
psd_mean = [];
figure()
for j = 2:NumChannels+1
    psd = [];
    for i = 1:size(StartPositions)
        psd(:,i) = pwelch(s(StartPositions(i):StopPositions(i),j),32,16,[],512);  
    end
    length(psd)
    psd_mean(j,:) = mean(psd,2);
    plot(log(psd_mean(j,:)))
    hold on
end


% For topoplots: we represent the PSD for specific time lapse: 
% What we can do is take the average PSD for each type

% For visualization: we can extract the power by squaring the filtered signal
% Then we can do a moving average


%% Averaging results per class:

% Extraction of data:
EventId = 773; % right cue

StartPositions = h1.EVENT.POS(h1.EVENT.TYP == EventId);
MinDuration = min(h1.EVENT.DUR(h1.EVENT.TYP == EventId))
StopPositions = StartPositions + h1.EVENT.DUR(h1.EVENT.TYP == EventId)-1;

NumTrials = length(StartPositions); % we get 30 trials
NumChannels = size(s,2) - 2;
Epoch = zeros(MinDuration, NumChannels, NumTrials);

% Epoching
for trial_id = 1:NumTrials
    cstart = StartPositions(trial_id);
    cstop = cstart + size(Epoch,1)-1;
    disp(['Right cue for trial ', num2str(trial_id), ' start at ', num2str(cstart), ' and stop at ', num2str(cstop)])
    Epoch(:,:,trial_id) = s(cstart:cstop, 2:NumChannels+1);
end

