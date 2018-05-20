clc
clear all
close all

% Which files do you want to analyse ?
files = 'Elise'; % Anonymous, Elise, Mathieu, Arnaud
% Do you want to compute the psd from raw data ?
create_psd = false;

% Choose the type of figure to plot: 
%   1 -> plot
%   0 -> no plot
topoplot = 0;
PSD_time_plot = 1;
Frequency_time_plot = 1;
Discrimancy_map = 1;

%% Create PSD matrices:

if create_psd 
    % Parameters:
    beta_band=0; %1=on,0=off
    mu_band =0;
    CAR=0;
    SmallLaplacian=0;
    LargeLaplacian=0;
    DoLabel=1;
    % do labe==1 if you want to do psd

    GetPSD(beta_band, mu_band, CAR, SmallLaplacian, LargeLaplacian, DoLabel, files)
end

%% Initialization:

addpath(genpath('functions'))

frequencies = load(['SPD/' files '/Frequences.mat']);
window_label = load(['SPD/' files '/WindowLabel.mat']);

% Load PSD estimates
psd_small_laplacian = load(['SPD/' files '/SPD with SmallLaplacian Spatial filtre.mat']);
psd_large_laplacian = load(['SPD/' files '/SPD with LargeLaplacian Spatial filtre.mat']);
psd_CAR_filter = load(['SPD/' files '/SPD with CAR Spatial filtre.mat']);
psd_no_spatial_filter = load(['SPD/' files '/SPD with NO Spatial filtre']);
psd = {psd_small_laplacian, psd_large_laplacian, psd_CAR_filter, psd_no_spatial_filter}; %ALL PSD

% Define interesting frequency bands
mu_band = [3:6];
beta_band = [7:12];
mu_beta_band = [3:12];
band = {mu_band, beta_band};

% Define psd window frequency
window_frequency = 16;

% Load the desired run
run = 1;
load(['SPD/' files '/WindowLabelRun.mat']);
if run > size(WindowLabelRun,1)
    error('The run you have chosen does not exist')
else
    run_ind = [WindowLabelRun(run,1):WindowLabelRun(run,2)];
end

load(['SPD/' files '/Event Window.mat'])
% Take only windows that correspond to the selected run
Event_window = Event_window(find((Event_window(:,2) > run_ind(1)) & (Event_window(:,5) < run_ind(end))),:);

% Choose the type of psd and frequency band 
%psd_selected = psd_large_laplacian.psdt(run_ind,:,:);
psd_selected = psd_large_laplacian.psdt;
band_selected = mu_beta_band;
name = 'Large Laplacian';

 

%% Get Epoching
[Epoch_both_feet, Epoch_both_hands, Baseline_both_feet, Baseline_both_hands, trial_length_feet, trial_length_hand] = Epoching(psd_selected, band_selected, Event_window);

%% Get topoplots

if topoplot
    all_action = false;
    all_plots = true;
    if all_plots
        psd_selected = psd;
    end
    GetTopoplot(psd_selected, band_selected, window_frequency, frequencies, Event_window, window_label, all_plots, all_action, name)
end

%% PSD time plot

if PSD_time_plot
    GetPSD_TimePlot(psd_selected, window_label.labelAction, frequencies.Frequencies, band_selected, window_label, files, Event_window);
end

%% Frequency time plots:

if Frequency_time_plot
    GetFrequencyTimePlot(psd_selected, band_selected, window_frequency, frequencies, files, Event_window);
end

%% Discrimancy maps:

if Discrimancy_map
    discrimancy = GetDiscrimancyMap(psd_selected, band_selected, window_frequency, frequencies, files, Event_window);
end
