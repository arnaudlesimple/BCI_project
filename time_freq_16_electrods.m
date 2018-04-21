clc
clear all
close all
% faire p-e une fonction
%% Extract data
load('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/SPD with CAR Spatial filtre.mat');
load('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/WindowLabel.mat');
load ('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/Frequences.mat');
load ('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/Feet Start_End window.mat');

%parameters
mu_band=1;
beta_band=0;
%if one is 0, the other must be 1

if mu_band
   band = [3:6];
end

if beta_band
   band=[7:18]; 
end

gap_percentage = 0.1;
% PSD_both_feet = mean(elec1(Start_End_Feet(1,1):Start_End_Feet(1,2), Frequencies(mu_band), :),2);%mean sur les frequences
% 10percent = (Start_End_Feet(1,2)-Start_End_Feet(1,1))
% PSD_BF_time = elec1(Start_End_Feet(1,1)-:,1,1);
% 
% index = psdt(labelAction==771, Frequencies(mu_band), :)
% PSD_both_feet_10percent = psdt(labelAction==771, Frequencies(mu_band), :)
% PSD_both_hand = mean(mean(psdt(labelAction==773, Frequencies(mu_band), :),1),2);


%% ARNAUD
load('SPD/Feet Start_End window.mat');
load('SPD/Hand Start_End window.mat');

num_trial = length(Start_End_Feet);
trial_length = min(Start_End_Feet(:,2)-Start_End_Feet(:,1)); %feet &hand length is the same


%have "rest values" before the trial
% zero_gap = gap_percentage*trial_length;
% zero_gap_vec = zeros(zero_gap)

%initialization
number_electrode = 16;
PSD_both_feet = zeros(trial_length,length(band),num_trial);

for n_electrode = 1:number_electrode
    for trial_number = 1:num_trial
        feet_trial_samples = [Start_End_Feet(trial_number,1):Start_End_Feet(trial_number,1) + trial_length-1];
        hand_trial_samples = [Start_End_Hand(trial_number,1):Start_End_Hand(trial_number,1) + trial_length-1];

        Epoch_both_feet(:,:,trial_number, n_electrode) = psdt(feet_trial_samples, band, n_electrode);
        % Extract mean PSD over time from both hand
        Epoch_both_hands(:,:,trial_number, n_electrode) = psdt(hand_trial_samples, band, n_electrode);  
    end
end
%% 
% electrode 10:
PSD_both_feet = mean(mean(Epoch_both_feet(:,:,:,10),3),2);
PSD_both_hands = mean(mean(Epoch_both_hands(:,:,:,10),3),2);

figure()
title('Electrode 10')
plot(log10(PSD_both_feet)); hold on;
plot(log10(PSD_both_hands));
hold on;
legend('both feet','both hands');


%% Plotting data

% en dB : 10log(freq)








