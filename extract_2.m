clc
clear all
close all
%%
% Hi all,
% I would like to provide some general steps for the analysis that you have to do.
% On Thursday we can discuss them.
% Grand average analysis
% % The goalis to highlight the neural correlates of the mental task performed by the subject(s) at the grand average level.
% In our case, the classical neural correlates for Motor Imagery tasks are the ERD/ERS (event related de-synchronization/syncronization)
% of the channels over the (pre-)motor cortex in the mu/beta bands.
% In the slides, you have the formula to compute the ERD/ERS and, as you know, this is based on the power of the signals during the activity period
% normalized with the power in the baseline period. In our protocol we can use the fixation as baseline (reference) period.
% To compute the power of the signals, you can use different approaches. The classical one is to compute the Power Spectral Density (PSD).
% To visualize the results you have different options: show the time series of the power in given frequency bands for all channels, topoplots by selecting
% specific window of the activity period (or even the whole period), using spectro-temporal visualization for each channel.
% Again the goal is to show that there are differences in ERD/ERS at the grand average level between the two motor imagery tasks performed by the subject(s).
% The rationale is that before starting to create a classifier you need to be sure that you can see differences at least at the grand average level.
% 
% 
% 
% 
% 4) Save the processed data for each file
% The processing is done, so it is better to save it (so you don't need to recompute every time the PSD).
% 5) Import and concatenate each processed file
% You can import the processed data of each file and than concatenate it. As you can see, the structure of the resulting data 
% is not dependent by the number of files you used. 
% Remember also to properly concatenate the position of the events.
% Dimension of the data: [NTotalWindows x NFrequencies x NChannels]
% 6) Create the vectors of labels
% Now, it's time to label your data. The idea is that you can create simple vectors that labels the first dimension of your data matrix (windows). 
% You can create different kind of these vectors.
% For sure you will need one that indicates the activity period and the class performed in this period. The vector labels will have the same length 
% of the first dimension of your data matrix (e.g., NTotalWindows). You need to use the events position, duration and type in order to create these
% vectors. Let's imagine that you create a vector label Ck corresponding to the activity period. 
% This vector can have 3 different values: 0-> no activity period, 771->activity period for both feet, 773->activity period for both hands. 
% If you need to take all the data related to class 771, then you need just to use this vector: PSD_both_feet = data(Ck==771, :, :).
% 
% 7) Extract the period of interest from the data
% By exploiting the vectors of labels, you can extract the data that you are interested in. Of course, 
% this depends on what you want to visualize. For instance, if you want to have the grand average PSD over the whole channels for the two classes, 
%     you do as follows:
% PSD_both_feet = squeeze(mean(data(Ck==771, :, :))):
% PSD_both_hands = squeeze(mean(data(Ck==773, :, :)));
% * squeeze is used to remove the first singleton dimension after the mean
% 
% 8) Visualize the data
% According to what you want to visualize, you can manipolate the data. For instance, for the topoplot you can show only the values of the channels (no temporal or spectral information). In this case, you can decide to average the data across the activity period and in the MuBand:
% topoplot(squeeze(squeeze(mean( mean(data(Ck==771, MuFreqIds, :), 1), 2 ) )))
% topoplot(squeeze(squeeze(mean( mean(data(Ck==773, MuFreqIds, :), 1), 2 ) )))
% Note that this is just possible visualization.
% 
% It is worth to notice that the processing flow can be different. I suggested you to use this one because can be easily re-used in the case of the Single-Trial analysis. In particular, you can use the same code until step 7. In any case, we will talk about single trial analysis in the next classes.
% Have a nice weekend,
% Luca

%% Parameter

beta_band=0; %1=on,0=off
mu_band =0;
CAR=0;
SmallLaplacian=1;
LargeLaplacian=0;

if CAR+SmallLaplacian+LargeLaplacian >1 | beta_band+mu_band >1
    A='Error'
    return;
end
%% Data extraction:
%1) For each gdf file, import data and events
%Dimension of the data: [NSamples x NChannels]

addpath(genpath('projects_common_material/biosig'));
filename1 = 'project2-data-example/anonymous.20170613.161402.offline.mi.mi_bhbf.gdf';
filename2 = 'project2-data-example/anonymous.20170613.162331.offline.mi.mi_bhbf.gdf';
filename3 = 'project2-data-example/anonymous.20170613.162934.offline.mi.mi_bhbf.gdf';
[rec1, info1] = sload(filename1);
[rec2, info2] = sload(filename2);
[rec3, info3] = sload(filename3);

% put data from 3 files into one matrix
AllRec = [rec1;rec2;rec3];
sizeRec = [length(rec1),length(rec2),length(rec3)];
Reclabel = [ones(sizeRec(1),1); ones(sizeRec(2),1)*2; ones(sizeRec(3),1)*3];

AllRecInit = [Reclabel, AllRec(:,[1:16])]; %without reference
% s: rows correspond to the samples
%    first column is the filename label
%    columns 2 to 17 -> 16 channels
%    column 18 -> reference electrode %je l'enl�ve

% Update position data
info2.EVENT.POS = info2.EVENT.POS + length(rec1);
info3.EVENT.POS = info3.EVENT.POS + length(rec2) + length(rec1);

%% Isolate interesting bandpass


% isolate mu and beta bandpass
sample_rate = info1.SampleRate;

if mu_band ==1
% mu band = [7-14] Hz
[b,a] = butter(4,[7,14]/sample_rate/2);
fvtool(b,a)
BandData = filter(b,a,AllRec(:,2:end)); %on prend pas colonne avec label
AllRecBand=[Reclabel,BandData];
end

if beta_band == 1
% beta band = [15-35] Hz
[b,a] = butter(4,[15,35]/sample_rate/2);
fvtool(b,a)
BandData = filter(b,a,AllRec(:,2:end));
AllRecBand=[Reclabel,BandData];
end

if beta_band ==0 && mu_band==0
    AllRecBand=AllRecInit;
end

%% Spatial filtering:

%2) For each gdf file, apply the spatial filter (laplacian or CAR)
%The spatial filter we have seen in class are time-independent. It means that they can be applied sample by sample (you don't need the past or the future). 
%Dimension of the data: [NSamples x NChannels]


% CAR filtering:
%chaque filtre prend AllRec contenant [label, data], freq peut deja �tre
%filtr�

if CAR==1   
    car_mu_data = [];
    for i = 1:length(AllRecBand)
        moy = mean(AllRecBand(i,2:end)); %moyenne sur tout les canaux � chaque T
        for j = 2:size(AllRecBand,2)
            car_mu_data(i,j) = AllRecBand(i,j) - moy;
            
        end
    end
    FilteredData=[Reclabel,car_mu_data];
end

% Small Laplacian filtering
if SmallLaplacian ==1   
    small_laplacian = load('small_laplacian.mat');
    FilteredData = [Reclabel, AllRecBand(:,2:end) * small_laplacian.lap];
end

% Large laplacian filter
if LargeLaplacian ==1
    large_laplacian = load('large_laplacian.mat');
    FilteredData = [Reclabel , AllRecBand(:,2:end) * large_laplacian.lap];
end

%filteredDate= [label,datafiltr�+band]
%% Power spectral density 

% 3) For each gdf file, compute the neural correlates of interest
% In the classical case where we want to extract the PSD, we can apply the pwelch (there are also other functions in matlab that 
% you can explore, for instance spectrogram).
% If you apply this kind of functions, you need a temporal window (-> you cannot compute the PSD on a single sample).
% In order to do that you need to window your signal. The easiest way
% is to create two vectors with the starting and stopping positions of all the window. Usually the window's length is 1s and the overlap 
% between windows is 0.0625 ms. 
% Then, for each window you compute the PSD (let's limit the frequency range from 4 Hz to 48 Hz, with frequency resolution of 2 Hz). 
% Moreover, I suggest you to compute it over the whole files (not only on the interest periods). 
% In this step, you are manipulating the temporal dimension (actually, you are downsampling 
% your signal [from samples->to window]) and you are adding a new dimension to your data in the spectral domain (NFrequencies). 
% On Thursday we will analyses more in detail this step.
% Dimension of the data: [NWindows x NFrequencies x NChannels]
% Remember, that you have also to convert the event positions and durations from samples to windows

% First, we estimate the PSD on the whole recording.
% We use windows of 1s, with a time overlap of 0.9375 s
% Thus we will get psd estimation of 1s, with the time varying of 62.5 ms at each estimation

time_window = 1; % in seconds
Init_freq = 512;
window_change=0.0625;
time_overlap = 1 - window_change; % in seconds


frequ_psdt=1/window_change;

% 62.5 ms = 16 samples
Frequencies=[4:2:48];

temps=size(FilteredData,1)/Init_freq;
sample=temps*frequ_psdt-(frequ_psdt-1); %car pour les dernieres on ne peut plus d�caler!!
psdt1=zeros(length(Frequencies),sample,16);
psdt=zeros(sample,length(Frequencies),16);
%%
window_start=zeros(sample,1);
window_end=zeros(sample,1);

for i=1:sample
    
window_start(i)=1+(i-1)*(window_change * Init_freq);% = depart absolu*#window*deltaSample
window_end(i)=window_start(i)+(time_window*Init_freq)-1;
end
%%

for n_electrode = 1:16
    
  [psdt1(:,:,n_electrode), w, t] = spectrogram(FilteredData(:,n_electrode+1),time_window * Init_freq, time_overlap * Init_freq, Frequencies ,512);

  psdt(:,:,n_electrode)=transpose(psdt1(:,:,n_electrode));
end

% psdt: columns = frequencies (there are 257 frequencies)
%       rows = time: each column represents 1 s but we have to take into account the 62.5 ms overlap
%       -> psdt(:,1) corresponds to PSD estimation between 0 s and 1 s
%       -> psdt(:,2) corresponds to PSD estimation between 0.0625 s and 1.0625 s
%       -> psdt(:,3) corresponds to PSD estimation between 0.1250 s and 1.1250 s
% That is why the number of columns is equal to
% number_of_samples_of_the_whole_signal * frequ_window/frequ_init - (frequ_window-1)

%% Epoching
% 4) Save the processed data for each file
% The processing is done, so it is better to save it (so you don't need to recompute every time the PSD).
% 5) Import and concatenate each processed file
% You can import the processed data of each file and than concatenate it. As you can see, the structure of the resulting data 
% is not dependent by the number of files you used. 
% Remember also to properly concatenate the position of the events.
% Dimension of the data: [NTotalWindows x NFrequencies x NChannels]

% 6) Create the vectors of labels
% Now, it's time to label your data. The idea is that you can create simple vectors that labels the first dimension of your data matrix 
% (windows). 
% You can create different kind of these vectors.
% For sure you will need one that indicates the activity period and the class performed in this period. 
% The vector labels will have the same length 
% of the first dimension of your data matrix (e.g., NTotalWindows). You need to use the events position, duration and type 
% in order to create these
% vectors. Let's imagine that you create a vector label Ck corresponding to the activity period. 
% This vector can have 3 different values: 0-> no activity period, 771->activity period for both feet, 773->activity period for 
% both hands. 
% If you need to take all the data related to class 771, then you need just to use this vector: PSD_both_feet = data(Ck==771, :, :).
% 


%Fixation = 786
%Both Hands/Right Cue = 773;
%Both Feet/Left Cue = 771;
%Continuous Feedback = 781;
%Boom Target Hit = 897;
%Boom Target Miss = 898;

% Extraction of data:
EventId = 773; % right cue

StartPositions = (info1.EVENT.POS(info1.EVENT.TYP == EventId));
EndPositions = StartPositions+(info1.EVENT.DUR(info1.EVENT.TYP == EventId));
Window_event773=zeros(1,length(StartPositions));

for i=1:length(StartPositions)
Debut(i)=find(window_start > StartPositions(i),1,'first')-1; %donne index donc deja le numero de la window!!
Fin(i)=find(window_end > EndPositions(i),1,'first')-1;
Window_event773(i)=[Debut(i):Fin(i)];

end

%%



NumTrials = length(StartPositions); 
NumFrequencies = size(psdt,1);
Epoch = zeros(NumFrequencies, EndPositions, NumTrials);

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

