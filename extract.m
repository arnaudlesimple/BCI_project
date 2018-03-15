clc
clear all
close all
addpath(genpath('biosig'));
filename1 = 'project2-data-example\anonymous.20170613.161402.offline.mi.mi_bhbf.gdf';
filename2 = 'project2-data-example\anonymous.20170613.162331.offline.mi.mi_bhbf.gdf';
filename3 = 'project2-data-example\anonymous.20170613.162934.offline.mi.mi_bhbf.gdf';
[s1, h1] = sload(filename1);
[s2, h2] = sload(filename2);
[s3, h3] = sload(filename3);

s = [s1;s2;s3];
size_ = [length(s1),length(s2),length(s3)];
label = [ones(size_(1),1); ones(size_(2),1)*2; ones(size_(3),1)*3];
s = [label, s];

h2.EVENT.POS = h2.EVENT.POS + length(s1);
h3.EVENT.POS = h3.EVENT.POS + length(s2) + length(s1);

%% Isolate interesting bandpass

%faire un filtre et isole les bandes mu et beta
sample_rate = h1.SampleRate;

% Bande mu = [7-14] Hz
[b,a] = butter(4,[7,14]/sample_rate/2);
fvtool(b,a)
mu_data = filter(b,a,s(:,1:end-1));

% Bande beta = [15-35] Hz
[b,a] = butter(4,[15,35]/sample_rate/2);
fvtool(b,a)
beta_data = filter(b,a,s);

%% Spatial filtering:

% CAR filtering:
car_mu_data = [];
for i = 1:length(mu_data)
    moy = mean(mu_data(i,:));
    for j = 1:size(mu_data,2)
        car_mu_data(i,j) = mu_data(i,j) - moy;
    end
end
%%
% Small Laplacian filtering
map_channels_small = {[4],[3,7],[2,4,8],[1,3,5,9],[4,6,10],[5,11],[2,8,12],[3,7,9,13],1,1,1,1,1,1,1,1};

for i = 1:length(mu_data)
    for j = 1:(size(mu_data,2)-1)
        moy = mean(mu_data(i,map_channels_small{j}));
        slap_mu_data(i,j) = mu_data(i,j) - moy;
    end
end

%%
% Large Laplacian filtering
map_channels_large = {};

for i = 1:length(mu_data)
    for j = 1:(size(mu_data,2)-1)
        moy = mean(mu_data(i,map_channels_large{j}));
        llap_mu_data(i,j) = mu_data(i,j) - moy;
    end
end