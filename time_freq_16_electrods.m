clc
clear all
close all
% faire p-e une fonction
%% Extract data
load('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/SPD with CAR Spatial filtre.mat');
load('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/WindowLabel.mat');
load ('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/Frequences.mat');

% struct pour 1 electrode electrode
elec1 = psdt(:, :, 1);
mu_band=[3:6];

PSD_both_feet = mean(elec1(:, Frequencies(mu_band), :),2); %mean sur les frequences
index_start_event = find (labelAction==both_feet,1,'first')

index = psdt(labelAction==771, Frequencies(mu_band), :)
PSD_both_feet_10percent = psdt(labelAction==771, Frequencies(mu_band), :)
PSD_both_hand = mean(mean(psdt(labelAction==773, Frequencies(mu_band), :),1),2);



%% Frequencies VS Baseline

% rapport frequency/baseline

%%
% mettre xlabel, ylabel, nom de la scale bar à ERD / ERS


%% Plotting data

% en dB : 10log(freq)








