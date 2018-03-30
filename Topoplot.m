clc
clear all
close all

addpath(genpath('projects_common_material/eeglab_current'))
addpath(genpath('projects_common_material/biosig'))
Map=load('projects_common_material/channel_location_16_10-20_mi.mat')
 filename1 = 'project2-data-example/anonymous.20170613.161402.offline.mi.mi_bhbf.gdf';
% filename2 = 'project2-data-example/anonymous.20170613.162331.offline.mi.mi_bhbf.gdf';
% filename3 = 'project2-data-example/anonymous.20170613.162934.offline.mi.mi_bhbf.gdf';
 [s1, h1] = sload(filename1);
% [s2, h2] = sload(filename2);
% [s3, h3] = sload(filename3);


B=topoplot(s1(1000,:),Map.chanlocs16);
% 
% A=5