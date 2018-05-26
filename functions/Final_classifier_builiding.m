clear
close all


ActionM = load('SPD/Elise/Event Window.mat');
ActionM = ActionM.Event_window;
EventM=load('SPD/Elise/WindowLabel.mat');%WindowLabel
EventM=EventM.labelAction;
Run=load('SPD/Elise/WindowLabelRun.mat');%WindowLabel
Run=Run.WindowLabelRun;

ActionTrain=ActionM(1:60,:);
Action1=ActionM(1:30,:);
Action2=ActionM(31:60,:);

ActionTest=ActionM(61:end,:);

EventTrain=EventM(1:14193,1);
Event1=EventM(1:6722,1);
Event2=EventM(6723:13266,1);

EventTest=EventM(13267:end,1);


%% Attention alpha doit etre penser en terme de practicit� car sinon on a pas le temps
psd_small_laplacian = load('SPD/Elise/SPD with SmallLaplacian Spatial filtre.mat');
psd_large_laplacian = load('SPD/Elise/SPD with LargeLaplacian Spatial filtre.mat');
psd_CAR_filter = load('SPD/Elise/SPD with CAR Spatial filtre.mat');
psd_no_spatial_filter = load('SPD/Elise/SPD with NO Spatial filtre');

selected_data=psd_small_laplacian.psdt;

% selected_data1=selected_data(1:6722,:,:);
% selected_data2=selected_data(6723:13266,:,:);
% selected_data=selected_data(1:13266,:,:);
% selected_data_Train=selected_data(13267:end,:,:);


window_frequency = 16;
frequencies = load('SPD/Anonymous/Frequences.mat');
load('SPD/Anonymous/Frequences.mat');

mu_band = 3:6;
beta_band = 7:18;
mu_beta_band = 3:18;
all_band=1:23;

band = {mu_band, beta_band};
band_selected = all_band;

Classifier=[1,2,3,4] %["linear","diaglinear","quadratic","diagquadratic"];
%%
discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies,'inutile',Action1);
discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies,'inutile',Action2);
%%
discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies,'inutile',ActionTrain);



%% Mathieu
features= [22 13; 24 15; 32 13; 30 13; 8 13; 22 6; 30 7; 6 7]; %[frequ x channel]
%% ELise
features=  [12 7; 10 7; 12 11; 10 11;16 7; 18 9; 32 4];
%%
features=sortrows(features,2);

alpha=0.06;
limit771=0.8;
limit773=0.2;
mean_trial=zeros(4,1);

Feat_selected=zeros(4,length(features));
%% test 4 classifier
%% FFS features

for Classifier=1:4
     [Feat_selected(Classifier,:),historyCrit,historyIn] = FFS(selected_data,ActionTrain,features,alpha,limit771,limit773,Frequencies,Classifier);
     mean_trial(Classifier)=historyCrit(end);
end

%% d
Final= [mean_trial,Feat_selected]

FFS_features={features(find(Feat_selected(1,:)),:),features(find(Feat_selected(2,:)),:),features(find(Feat_selected(3,:)),:),features(find(Feat_selected(4,:)),:)};


%% different alpha

alpha=0.01:0.01:0.07;
mean_window=zeros(length(alpha),4);
mean_trial=zeros(length(alpha),4);



for Classifier=1:4


    for i=1:length(alpha)

        [mean_window(i,Classifier),std_window,mean_trial(i,Classifier),std_trial,mean_Confusion]=GetClassifierAccuracy(selected_data,ActionTrain,FFS_features{Classifier},alpha(i),limit771,limit773,Frequencies,Classifier);

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
xlabel('Alpha [ ]')
ylabel('Test error [%]')

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
xlabel('Alpha [ ]')
ylabel('Test error [%]')

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
xlabel('Alpha [ ]')
ylabel('Test error [%]')

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
xlabel('Alpha [ ]')
ylabel('Test error [%]')

%% confusion matrix of optimal of each optimal classifier

mean_Confusion=zeros(4,3,3);

for Classifier=1:4

        [mean_window(i,Classifier),std_window,mean_trial(i,Classifier),std_trial,mean_Confusion(Classifier,:,:)]=GetClassifierAccuracy(selected_data,ActionTrain,FFS_features{Classifier},OptimalAlpha(Classifier),limit771,limit773,Frequencies,Classifier);
 
end



%% test on new data
Param={OptimalAlpha;FFS_features}
Classifier=3;


[avg_s_s_error,STD_s_s_error,avg_trial_error,STD_trial_error,FINALClassifier]=ClassifierAccuracyOnNewSample(selected_data,ActionTest,ActionTrain,FFS_features{Classifier},OptimalAlpha(Classifier),limit771,limit773,Frequencies,Classifier)

%%
freq=zeros(size(Final_features,1),1);

for d=1:size(Final_features,1)
 freq(d)=find(Frequencies==Final_features(d,1));
 end

selected_features=[freq,Final_features(:,2)]

%%
small_laplacian = load('small_laplacian.mat');
small_laplacian=small_laplacian.lap;

param=[OptimalAlpha(Classifier),OptimalLim(Classifier)]

%%
save('Classifier\supportMathieu','FINALClassifier','small_laplacian','param','selected_features')






% function [class_error]=error(real_label, predicted_label)
%     false=nnz(real_label-predicted_label);
%     
%     class_error=(false/length(real_label));
% end

