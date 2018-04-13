clc
clear all
close all

%% Data extraction:

addpath(genpath('projects_common_material/biosig'));
filename1 = 'project2-data-example/anonymous.20170613.161402.offline.mi.mi_bhbf.gdf';
filename2 = 'project2-data-example/anonymous.20170613.162331.offline.mi.mi_bhbf.gdf';
filename3 = 'project2-data-example/anonymous.20170613.162934.offline.mi.mi_bhbf.gdf';
[s1, h1] = sload(filename1);
[s2, h2] = sload(filename2);
[s3, h3] = sload(filename3);

% put data from 3 files into one matrix
s = [s1;s2;s3];
size_ = [length(s1),length(s2),length(s3)];
label = [ones(size_(1),1); ones(size_(2),1)*2; ones(size_(3),1)*3];
s = [label, s];
% s: rows correspond to the samples
%    first column is the filename label
%    columns 2 to 17 -> 16 channels
%    column 18 -> reference electrode

% Update position data
h2.EVENT.POS = h2.EVENT.POS + length(s1);
h3.EVENT.POS = h3.EVENT.POS + length(s2) + length(s1);

%% Isolate interesting bandpass

% isolate mu and beta bandpass
sample_rate = h1.SampleRate;

% mu band = [7-14] Hz
[b,a] = butter(4,[7,14]/sample_rate/2);
fvtool(b,a)
mu_data = filter(b,a,s(:,2:end));

% beta band = [15-35] Hz
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
small_laplacian = load('small_laplacian.mat');
small_lap_mu_data = mu_data(:,1:end-1) * small_laplacian.lap;

% Large laplacian filter
large_laplacian = load('large_laplacian.mat');
large_lap_mu_data = mu_data(:,1:end-1) * small_laplacian.lap;

%% Power spectral density 

% First, we estimate the PSD on the whole recording.
% We use windows of 1s, with a time overlap of 0.9375 s
% Thus we will get psd estimation of 1s, with the time varying of 62.5 ms at each estimation

time_window = 1; % in seconds
freq_s = 512;
time_overlap = 1 - 0.0625; % in seconds
% 62.5 ms = 32 samples

for n_electrode = 1:16
    [psdt(:,:,n_electrode), w, t] = spectrogram(small_lap_mu_data(:,n_electrode),time_window * freq_s, time_overlap * freq_s, [7:14] ,512);
end

% psdt: rows = frequencies (there are 257 frequencies)
%       columns = time: each column represents 1 s but we have to take into account the 62.5 ms overlap
%       -> psdt(:,1) corresponds to PSD estimation between 0 s and 1 s
%       -> psdt(:,2) corresponds to PSD estimation between 0.0625 s and 1.0625 s
%       -> psdt(:,3) corresponds to PSD estimation between 0.1250 s and 1.1250 s
% That is why the number of columns is equal to
% number_of_samples_of_the_whole_signal / number_of_samples_of_62.5 ms

%% Epoching

%Fixation = 786
%Both Hands/Right Cue = 773;
%Both Feet/Left Cue = 771;
%Continuous Feedback = 781;
%Boom Target Hit = 897;
%Boom Target Miss = 898;

% Extraction of data:
EventId = 773; % right cue

StartPositions = round(h1.EVENT.POS(h1.EVENT.TYP == EventId) / 32);
MinDuration = round(min(h1.EVENT.DUR(h1.EVENT.TYP == EventId)) / 32);

NumTrials = length(StartPositions); 
NumFrequencies = size(psdt,1);
Epoch = zeros(NumFrequencies, MinDuration, NumTrials);

% Epoching
for n_electrode = 1:16
    for trial_id = 1:NumTrials
        cstart = StartPositions(trial_id);
        cstop = cstart + size(Epoch,2)-1;
        disp(['Right cue for trial ', num2str(trial_id), ' start at ', num2str(cstart), ' and stop at ', num2str(cstop)])
        Epoch(:,:,trial_id, n_electrode) = psdt(:,cstart:cstop, n_electrode);
    end
end

%% Averaging results 

% The spectrogram function returns complex numbers, which include both magnitude and phase information.
% we consider only the magnitude
figure()
for channel = 1:16
    AverageEpoch = mean(Epoch(:,:,:,channel),3);
    if channel <= 9
        subplot(3,3,channel)
    elseif channel == 10
        figure()    
        subplot(3,3,channel-9)
    else
        subplot(3,3,channel-9)
    end
    imagesc(abs(AverageEpoch))
    title(['right cue channel ' mat2str(channel)])
    ylabel('frequency')
    xlabel('time')
end


%%
figure()
W_unnorm=W*512;
plot(W_unnorm,log(psd))
% we see alpha peak
% We must check it looks like this for our own signals.


