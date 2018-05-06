clear
close all


Action = load('SPD/Event Window.mat');
Action = Action.Event_window;
Event=load('SPD/WindowLabel.mat');%WindowLabel
Event=Event.labelAction;

%%
psd_small_laplacian = load('SPD/SPD with SmallLaplacian Spatial filtre.mat');
psd_large_laplacian = load('SPD/SPD with LargeLaplacian Spatial filtre.mat');
psd_CAR_filter = load('SPD/SPD with CAR Spatial filtre.mat');
psd_no_spatial_filter = load('SPD/SPD with NO Spatial filtre');

selected_data=psd_CAR_filter.psdt;

window_frequency = 16;
frequencies = load('SPD/Frequences.mat');
load('SPD/Frequences.mat');

mu_band = 3:6;
beta_band = 7:18;
mu_beta_band = 3:18;
all_band=1:23;

band = {mu_band, beta_band};
band_selected = all_band;
%%
discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies);


%% different alpha

features= [24 9; 12 2; 12 3; 12 5; 12 7;12 8;12 11]; %[frequ x channel]

alpha=0:0.1:1;
limit771=0.7;
limit773=0.3;

mean_window=zeros(1,length(alpha));
mean_trial=zeros(1,length(alpha));

for i=1:length(alpha)
    
    [mean_window(i),std_window,mean_trial(i),std_trial]=GetClassifierAccuracy(selected_data,Action,features,alpha(i),limit771,limit773) 
end

plot(alpha,mean_window,'r')
hold on;
plot(alpha,mean_trial,'b')
    

%% different number of features

%discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies);
features= [32 13;4 3;24 9; 12 2; 12 3; 12 5; 12 7;12 8;12 11;12 6; 12 10; 10 9; 12 4; 16 2 ;10 5;4 15;4 3; 12 13; 14 11; 24 16;30 2;30 12]; %[frequ x channel]


alpha=0.2;
limit771=0.7;
limit773=0.3;

mean_window=zeros(1,length(features));
mean_trial=zeros(1,length(features));

for i=2:length(features)
    
    [mean_window(i),std_window,mean_trial(i),std_trial]=GetClassifierAccuracy(selected_data,Action,features(1:i,:),alpha,limit771,limit773) 
end

plot(1:length(features),mean_window,'r')
hold on;
plot(1:length(features),mean_trial,'b')
    
%% same but with PCA

%discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies);
features= [24 9; 12 2; 12 3; 12 5; 12 7;12 8;12 11;12 6; 12 10; 10 9; 12 4; 16 2 ;10 5;4 15]; %[frequ x channel]


alpha=0.2;
limit771=0.7;
limit773=0.3;

mean_window=zeros(1,length(features));
mean_trial=zeros(1,length(features));

for i=2:length(features)
    
    [mean_window(i),std_window,mean_trial(i),std_trial]=GetClassifierAccuracy(selected_data,Action,features(1:i,:),alpha,limit771,limit773) 
end

plot(1:length(features),mean_window,'r')
hold on
plot(1:length(features),mean_trial,'b')
    