clc
clear all
close all
%% Extract data
load('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/SPD with CAR Spatial filtre.mat');
load('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/WindowLabel.mat');
load ('/Users/kientzelise/Documents/EPFL/Master/Master2/BCI/Projet/ProjetGit/BCI_project/SPD/Frequences.mat');

tot_num_elec = size(psdt,3);
band = 'mu_band'; 


for n_electrode=1:tot_num_elec
    figure(1)
    subplot(4,4,n_electrode)
    plot_time_freq(psdt, labelAction, Frequencies, n_electrode, band); hold on;
    title(['Electrode ', num2str(n_electrode)])
end
legend('both feet','both hands','end of fixation');


