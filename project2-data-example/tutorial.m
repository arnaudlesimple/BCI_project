clc
clear all
close all
addpath(genpath('../biosig'));
filename = 'anonymous.20170613.161402.offline.mi.mi_bhbf.gdf';

[s, h] = sload(filename);
% s = eeg signal -> we have 16 signal
% h is the header that gives information of eeg data

h.EVENT;
% need to know when events are appearing : fixation, cue, continuous feedback, boom ...
% There is the type (reference number), the position, the duration (number of samples)
[h.EVENT.TYP, h.EVENT.POS, h.EVENT.DUR];
% h.EVENT.DUR => duration is in number of samples.

% Total time
size(s,1)/h.SampleRate;
%%
% Steps: Fixation, cue, and continuous feedback. We need to extract these
% events (start position, duration...)

% Take all the cue periods: codes are 771 and 773
start_pos_cue = h.EVENT.POS(h.EVENT.TYP == 771);
stop_pos_cue = start_pos_cue + h.EVENT.DUR(h.EVENT.TYP == 771)-1;

% 781 is the continuous feedback
start_pos_cont = h.EVENT.POS(h.EVENT.TYP == 781);
stop_pos_cont = start_pos_cont + h.EVENT.DUR(h.EVENT.TYP == 781)-1;

%% Extraction of data:
EventId = 781;

StartPositions = h.EVENT.POS(h.EVENT.TYP == EventId);
MinDuration = min(h.EVENT.DUR(h.EVENT.TYP == EventId))
StopPositions = StartPositions + h.EVENT.DUR(h.EVENT.TYP == EventId)-1;

NumTrials = length(StartPositions); % we get 30 trials
NumChannels = size(s,2) - 1;
Epoch = zeros(MinDuration, NumChannels, NumTrials);

% Epoching
for trial_id = 1:NumTrials
    cstart = StartPositions(trial_id);
    cstop = cstart + size(Epoch,1)-1;
    disp(['Continuous feedback for trial ', num2str(trial_id), ' start at ', num2str(cstart), ' and stop at ', num2str(cstop)])
    Epoch(:,:,trial_id) = s(cstart:cstop, 1:NumChannels);
end

%%
TrialLb = h.EVENT.TYP(h.EVENT.TYP == 771 | h.EVENT.TYP == 773);

AverageClass1 = mean(Epoch(:,9,TrialLb == 771),3);
AverageClass2 = mean(Epoch(:,9,TrialLb == 773),3);

figure()
subplot(1,2,1)
plot(AverageClass1)
title('Channel Cz (9): Class 771')
subplot(1,2,2)
plot(AverageClass2)
title('Channel Cz (9): Class 773')

%%
data = Epoch(:,9,1);
[b,a] = butter(4,[8,12]/512/2);
fvtool(b,a)
data_filtered = filter(b,a,data)
figure()
plot(data)
hold on
plot(data_filtered)

% Suggestion: filtering on the whole file, and then Epoch