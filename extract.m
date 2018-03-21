clc
clear all
close all
addpath(genpath('biosig'));
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
%%
% Small Laplacian filtering
laplacian = load('laplacian_16_10-20_mi.mat');
for i = 1:length(mu_data)
    for j = 1:(size(mu_data,2)-1) 
        sum_channel = 0;
        for z = 1:(size(laplacian.lap,2)-1) 
            sum_channel = sum_channel + laplacian.lap(j,z) * mu_data(i,z);
        end
        slap_mu_data(i,j) = sum_channel / sum(laplacian.lap(j,:));
    end
end

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

