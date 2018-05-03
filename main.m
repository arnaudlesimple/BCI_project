clc
clear all
close all

% Choose the analysis:
create_psd = false;

% Choose the type of figure to plot: 
%   1 -> plot
%   0 -> no plot
topoplot = 1;
PSD_time_plot = 0;
PSD_time_plot_Elise = 0;
Frequency_time_plot = 0;
Discrimancy_map = 0;



%% Create PSD matrices:

if create_psd 
    % Parameter
    beta_band=0; %1=on,0=off
    mu_band =1;
    CAR=1;
    SmallLaplacian=0;
    LargeLaplacian=0;
    DoLabel=0;

    GetPSD(beta_band, mu_band, CAR, SmallLaplacian, LargeLaplacian, DoLabel)
end

%% Initialization:

addpath(genpath('functions'))

frequencies = load('SPD/Frequences.mat');
window_label = load('SPD/WindowLabel.mat');
load('SPD/Event Window.mat')

% Load PSD estimates
psd_small_laplacian = load('SPD/SPD with SmallLaplacian Spatial filtre.mat');
psd_large_laplacian = load('SPD/SPD with LargeLaplacian Spatial filtre.mat');
psd_CAR_filter = load('SPD/SPD with CAR Spatial filtre.mat');
psd_no_spatial_filter = load('SPD/SPD with NO Spatial filtre');
psd = {psd_small_laplacian, psd_large_laplacian, psd_CAR_filter, psd_no_spatial_filter};

% Define interesting frequency bands
mu_band = [3:6];
beta_band = [7:18];
mu_beta_band = [3:18];
band = {mu_band, beta_band};

% Define psd window frequency
window_frequency = 16;


% Choose the type of psd and frequency band 
all_plots = false;
psd_selected = psd_large_laplacian;
band_selected = mu_beta_band;

%% Get Epoching
[Epoch_both_feet, Epoch_both_hands, Baseline_both_feet, Baseline_both_hands, trial_length_feet, trial_length_hand] = Epoching(psd_selected.psdt, band_selected);

%% Get topoplots

if topoplot
    all_action = true;
    GetTopoplot(psd_selected, band_selected, window_frequency, frequencies, Event_window, window_label, all_plots, all_action)
end

%% PSD time plot:

if PSD_time_plot
    GetPSD_TimePlotTry(Epoch_both_feet, Epoch_both_hands, Baseline_both_feet, Baseline_both_hands)
end

%% Frequency time plots:

if Frequency_time_plot
    GetFrequencyTimePlot(psd_selected.psdt, band_selected, window_frequency, frequencies);
end

%% Discrimancy maps:

if Discrimancy_map
    discrimancy = GetDiscrimancyMap(psd_selected.psdt, band_selected, window_frequency, frequencies);
end

%% PSD time plot: Elise w8

if PSD_time_plot_Elise
    tot_num_elec = size(psd_small_laplacian.psdt,3);
    band = 'mu_band'; 

    for n_electrode=1:tot_num_elec
        figure(1)
        subplot(4,4,n_electrode)
        [Epoch_both_feet, Epoch_both_hands] = GetPSD_TimePlot(psd_selected.psdt, window_label.labelAction, frequencies.Frequencies, n_electrode, band, window_label); hold on;
        title(['Electrode ', num2str(n_electrode)])
    end
    legend('both feet','both hands','end of fixation');
end