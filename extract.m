clc
clear all
close all
addpath(genpath('C:\Users\arnau\Documents\ETUDES\Master\Semestre_3\Brain Computer Interaction\Projects - Common material-20180301\biosig'));
filename1 = 'project2-data-example\anonymous.20170613.161402.offline.mi.mi_bhbf.gdf';
filename2 = 'project2-data-example\anonymous.20170613.162331.offline.mi.mi_bhbf.gdf';
filename3 = 'project2-data-example\anonymous.20170613.162934.offline.mi.mi_bhbf.gdf';
[s1, h1] = sload(filename1);
[s2, h2] = sload(filename2);
[s3, h3] = sload(filename3);

s = [s1;s2;s3];
size = [length(s1),length(s2),length(s3)];
label = [ones(size(1),1); ones(size(2),1)*2; ones(size(3),1)*3];
s = [label, s];

h2.EVENT.POS = h2.EVENT.POS + length(s1);
h3.EVENT.POS = h3.EVENT.POS + length(s2) + length(s1);



COUCOU
