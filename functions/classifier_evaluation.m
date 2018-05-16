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

Classifier=[1,2,3,4] %["linear","diaglinear","quadratic","diagquadratic"];
%%
discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies);
%%
features= [24 9; 12 2; 12 3; 12 5; 12 7;12 8;12 11]; %[frequ x channel]

alpha=0:0.1:1;
limit771=0.7;
limit773=0.3;


%% test 4 classifier
%% different alpha
mean_window=zeros(length(alpha),4);
mean_trial=zeros(length(alpha),4);

for Classifier=1:4


for i=1:length(alpha)
    
    [mean_window(i,Classifier),std_window,mean_trial(i,Classifier),std_trial]=GetClassifierAccuracy(selected_data,Action,features,alpha(i),limit771,limit773,Frequencies,Classifier);
end
end
OptimalAlpha=zeros(1,4);
%%
subplot(2,2,1)
plot(alpha,mean_window(:,1),'r')
hold on;
plot(alpha,mean_trial(:,1),'b')
[a,b]=min(mean_trial(:,1));
plot(alpha(b),a,'o','linewidth',1);
title('Linear classifier');
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal alpha =' num2str(alpha(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalAlpha(1)=alpha(b);

subplot(2,2,2)
plot(alpha,mean_window(:,2),'r')
hold on;
plot(alpha,mean_trial(:,2),'b')
[a,b]=min(mean_trial(:,2));
plot(alpha(b),a,'o','linewidth',1);
title('Diaglinear classifier')
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal alpha =' num2str(alpha(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalAlpha(2)=alpha(b);

subplot(2,2,3)
plot(alpha,mean_window(:,3),'r')
hold on;
plot(alpha,mean_trial(:,3),'b')
[a,b]=min(mean_trial(:,3));
plot(alpha(b),a,'o','linewidth',1);
title('Quadratic classifier')
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal alpha =' num2str(alpha(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalAlpha(3)=alpha(b);

subplot(2,2,4)
plot(alpha,mean_window(:,4),'r')
hold on;
plot(alpha,mean_trial(:,4),'b')
[a,b]=min(mean_trial(:,4));
plot(alpha(b),a,'o','linewidth',1);
title('Diagquadratic classifier')
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal alpha =' num2str(alpha(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalAlpha(4)=alpha(b);

%% different number of features

%discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies);
%features=[12 13; 14 11; 24 16; 12 10; 10 9; 12 4; 16 2 ;10 5;4 15;4 3]; %[frequ x channel]


alpha=0.2;
limit771=0.7;
limit773=0.3;

mean_window=zeros(length(features),4);
mean_trial=zeros(length(features),4);

for Classifier=1:4


for i=1:size(features,1)
    
   
    [mean_window(i,Classifier),std_window,mean_trial(i,Classifier),std_trial]=GetClassifierAccuracy(selected_data,Action,features(1:i,:),OptimalAlpha(Classifier),limit771,limit773,Frequencies,Classifier) ;
end
end
%%
close all;

xaxis=1:1:length(features);
OptimalFeature=zeros(1,4);


subplot(2,2,1)
plot(xaxis,mean_window(:,1),'r')
hold on;
plot(xaxis,mean_trial(:,1),'b')
[a,b]=min(mean_trial(:,1));
plot(xaxis(b),a,'o','linewidth',1);
title('Linear classifier');
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal #features =' num2str(xaxis(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalFeature(1)=xaxis(b);

subplot(2,2,2)
plot(xaxis,mean_window(:,2),'r')
hold on;
plot(xaxis,mean_trial(:,2),'b')
[a,b]=min(mean_trial(:,2));
plot(xaxis(b),a,'o','linewidth',1);
title('Diaglinear classifier')
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal #features =' num2str(xaxis(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalFeature(2)=xaxis(b);

subplot(2,2,3)
plot(xaxis,mean_window(:,3),'r')
hold on;
plot(xaxis,mean_trial(:,3),'b')
[a,b]=min(mean_trial(:,3));
plot(xaxis(b),a,'o','linewidth',1);
title('Quadratic classifier')
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal #features =' num2str(xaxis(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalFeature(3)=xaxis(b);

subplot(2,2,4)
plot(xaxis,mean_window(:,4),'r')
hold on;
plot(xaxis,mean_trial(:,4),'b')
[a,b]=min(mean_trial(:,4));
plot(xaxis(b),a,'o','linewidth',1);
title('Diagquadratic classifier')
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal #features =' num2str(xaxis(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalFeature(4)=xaxis(b);

%% test different limit


limit=0:0.05:0.5;
limit771=0.7;
limit773=0.3;

mean_window=zeros(size(limit,2),4);
mean_trial=zeros(size(limit,2),4);

for Classifier=1:4


for i=1:size(limit,2)
    
   
    [mean_window(i,Classifier),std_window,mean_trial(i,Classifier),std_trial]=GetClassifierAccuracy(selected_data,Action,features(1:OptimalFeature(Classifier),:),OptimalAlpha(Classifier),1-limit(i),limit(i),Frequencies,Classifier) ;
end
end
%%
close all;


OptimalLim=zeros(1,4);


subplot(2,2,1)
plot(limit,mean_window(:,1),'r')
hold on;
plot(limit,mean_trial(:,1),'b')
[a,b]=min(mean_trial(:,1));
plot(limit(b),a,'o','linewidth',1);
title('Linear classifier');
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal limit =' num2str(limit(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalLim(1)=limit(b);

subplot(2,2,2)
plot(limit,mean_window(:,2),'r')
hold on;
plot(limit,mean_trial(:,2),'b')
[a,b]=min(mean_trial(:,2));
plot(limit(b),a,'o','linewidth',1);
title('Linear classifier');
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal limit =' num2str(limit(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalLim(2)=limit(b);

subplot(2,2,3)
plot(limit,mean_window(:,3),'r')
hold on;
plot(limit,mean_trial(:,3),'b')
[a,b]=min(mean_trial(:,3));
plot(limit(b),a,'o','linewidth',1);
title('Linear classifier');
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal limit =' num2str(limit(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalLim(3)=limit(b);

subplot(2,2,4)
plot(limit,mean_window(:,4),'r')
hold on;
plot(limit,mean_trial(:,4),'b')
[a,b]=min(mean_trial(:,4));
plot(limit(b),a,'o','linewidth',1);
title('Linear classifier');
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal limit =' num2str(limit(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalLim(4)=limit(b);
