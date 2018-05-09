clc
clear all
close all

frequencies = load('SPD/Frequences.mat');
window_label = load('SPD/WindowLabel.mat');
load('SPD/Event Window.mat')

% Load PSD estimates
psd_small_laplacian = load('SPD/SPD with SmallLaplacian Spatial filtre.mat');
psd_large_laplacian = load('SPD/SPD with LargeLaplacian Spatial filtre.mat');
psd_CAR_filter = load('SPD/SPD with CAR Spatial filtre.mat');
psd_no_spatial_filter = load('SPD/SPD with NO Spatial filtre');
psd = {psd_small_laplacian, psd_large_laplacian, psd_CAR_filter, psd_no_spatial_filter};

mu_band = [3:6];
beta_band = [7:18];
band = {mu_band, beta_band};

window_frequency = 16;

%% Feature plots:

% obtain all plots:
for psd_num = 1:length(psd)
    for band_num = 1:length(band)
        GetFeaturePlot(psd{psd_num}, band{band_num}, frequencies, window_label, psd_num, band_num);
    end
end

% obtain 1 plot:
%GetFeaturePlot(psd_large_laplacian, mu_band, frequencies, window_label, psd_num, band_num);

%% Frequency time plots:

% For all plots:
for psd_num = 1:length(psd)
    for band_num = 1:length(band)
        GetFrequencyTimePlot(psd{psd_num}, band{band_num}, window_frequency, frequencies)
    end
end

% For 1 plot:
%GetFrequencyTimePlot(psd_small_laplacian, mu_band, window_frequency, frequencies)