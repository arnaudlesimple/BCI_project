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


%% Attention alpha doit etre penser en terme de practicité car sinon on a pas le temps
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
features=  [12 7; 10 7; 12 11; 10 11];

features=sortrows(features,2);

alpha=0.01:0.01:0.07;
limit771=0.8;
limit773=0.2;


%% test 4 classifier
%% different alpha
mean_window=zeros(length(alpha),4);
mean_trial=zeros(length(alpha),4);

for Classifier=1:4


for i=1:length(alpha)
    
    [mean_window(i,Classifier),std_window,mean_trial(i,Classifier),std_trial]=GetClassifierAccuracy(selected_data,ActionTrain,features,alpha(i),limit771,limit773,Frequencies,Classifier);
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

%% different number of features

%discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies);
%features=[12 13; 14 11; 24 16; 12 10; 10 9; 12 4; 16 2 ;10 5;4 15;4 3]; %[frequ x channel]



limit771=0.7;
limit773=0.3;

mean_window=zeros(length(features),4);
mean_trial=zeros(length(features),4);

for Classifier=1:4


for i=1:size(features,1)
    
   
    [mean_window(i,Classifier),std_window,mean_trial(i,Classifier),std_trial]=GetClassifierAccuracy(selected_data,ActionTrain,features(1:i,:),OptimalAlpha(Classifier),limit771,limit773,Frequencies,Classifier) ;
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
xlabel('# of features [ ]')
ylabel('Test error [%]')

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
xlabel('# of features [ ]')
ylabel('Test error [%]')

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
xlabel('# of features [ ]')
ylabel('Test error [%]')

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
xlabel('# of features [ ]')
ylabel('Test error [%]')


%% FFS features

Feat_selected=zeros(4,size(features,1));
mean_trial=zeros(4,1);

for Classifier=1:4
     [Feat_selected(Classifier,:),historyCrit,historyIn] = FFS(selected_data,ActionTrain,features,OptimalAlpha(Classifier),limit771,limit773,Frequencies,Classifier);
     mean_trial(Classifier)=historyCrit(end);
end

%% d
Final= [mean_trial,Feat_selected]

FFS_features={features(find(Feat_selected(1,:)),:),features(find(Feat_selected(2,:)),:),features(find(Feat_selected(3,:)),:),features(find(Feat_selected(4,:)),:)};

%% test different limit

limit=0:0.05:0.5;
limit771=0.7;
limit773=0.3;

mean_window=zeros(size(limit,2),4);
mean_trial=zeros(size(limit,2),4);

for Classifier=1:4


for i=1:size(limit,2)
    
   
    [mean_window(i,Classifier),std_window,mean_trial(i,Classifier),std_trial]=GetClassifierAccuracy(selected_data,ActionTrain,FFS_features{Classifier},OptimalAlpha(Classifier),1-limit(i),limit(i),Frequencies,Classifier) ;
end
end

%features(1:OptimalFeature(Classifier),:)
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
xlabel('threshold [%]')
ylabel('Test error [%]')

subplot(2,2,2)
plot(limit,mean_window(:,2),'r')
hold on;
plot(limit,mean_trial(:,2),'b')
[a,b]=min(mean_trial(:,2));
plot(limit(b),a,'o','linewidth',1);
title('Diaglinear classifier');
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal limit =' num2str(limit(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalLim(2)=limit(b);
xlabel('threshold [%]')
ylabel('Test error [%]')

subplot(2,2,3)
plot(limit,mean_window(:,3),'r')
hold on;
plot(limit,mean_trial(:,3),'b')
[a,b]=min(mean_trial(:,3));
plot(limit(b),a,'o','linewidth',1);
title('Quadratic classifier');
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal limit =' num2str(limit(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalLim(3)=limit(b);
xlabel('threshold [%]')
ylabel('Test error [%]')

subplot(2,2,4)
plot(limit,mean_window(:,4),'r')
hold on;
plot(limit,mean_trial(:,4),'b')
[a,b]=min(mean_trial(:,4));
plot(limit(b),a,'o','linewidth',1);
title('Diagquadratic classifier');
ylim([0,1]);
optimal=['minimal error is ' num2str(a) ' and optimal limit =' num2str(limit(b))];
legend('mean window class error','mean trial class error',optimal)
OptimalLim(4)=limit(b);
xlabel('threshold [%]')
ylabel('Test error [%]')


%% test on new data
Param=[OptimalAlpha;OptimalFeature;OptimalLim]
Classifier=3;

Final_features=features(1:OptimalFeature(Classifier),:);


[avg_s_s_error,STD_s_s_error,avg_trial_error,STD_trial_error,FINALClassifier]=ClassifierAccuracyOnNewSample(selected_data,ActionTest,ActionTrain,Final_features,OptimalAlpha(Classifier),1-OptimalLim(Classifier),OptimalLim(Classifier),Frequencies,Classifier)

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




%% 3D

limit=0:0.05:0.5;
alpha=0:0.1:1;

limit771=0.7;
limit773=0.3;

mean_window=zeros(4,size(limit,2),size(alpha,2));
mean_trial=zeros(4,size(limit,2),size(alpha,2));

for Classifier=1:1


    for i=1:size(limit,2)
           for j=1:size(alpha,2)
    
            [mean_window(Classifier,i,j),std_window,mean_trial(Classifier,i,j),std_trial]=GetClassifierAccuracy(selected_data,ActionTrain,features,alpha(j),1-limit(i),limit(i),Frequencies,Classifier) ;
           end
   end
end
%%

%subplot(2,2,1)
A=ones(11,11)-squeeze(mean_trial(1,:,:));

title('Linear classifier mean trial error');
xlabel('Limit');
ylabel('alpha');
zlabel('Accuracy');
t=squeeze(mean_trial(1,:,:));

[a1,b1]=min(t);
[a2,b2]=min(min(t));
b=b1(b2); %optimal_limit=b et b2=optimal alpha

%surf([limit(b)-0.025,limit(b),limit(b)+0.025],[alpha(b2)-0.05,alpha(b2),alpha(b2)+0.05],ones(3,3)*A(b,b2));

%hold on;
surf(limit,alpha,A);
optimal=['minimal error is ' num2str(a2) ' and optimal limit =' num2str(limit(b)) 'with alpha= ' num2str(alpha(b2))];
legend('mean trial class error',optimal)


%%
subplot(2,2,2)
A=ones(11,11)-squeeze(mean_trial(2,:,:));
surf(limit,alpha,A);
title('Diaginear classifier mean trial error');
xlabel('Limit');
ylabel('alpha');
zlabel('Accuracy');

subplot(2,2,3)
A=ones(11,11)-squeeze(mean_trial(3,:,:));
surf(limit,alpha,A);
title('Quadratic classifier mean trial error');
xlabel('Limit');
ylabel('alpha');
zlabel('Accuracy');

subplot(2,2,4)
A=ones(11,11)-squeeze(mean_trial(4,:,:));
surf(limit,alpha,A);
title('Diagquadratic classifier mean trial error');
xlabel('Limit');
ylabel('alpha');
zlabel('Accuracy');



Final_features_index=find(Feat_selected(a,:));
Final_features=features(Final_features_index,:)

%% build final classifier
    a=1;
     
    
    %%
    %Je recupere les temps de chaque trial ainsi que son label
start_feedback=ActionTrain(:,4);
end_feedback=ActionTrain(:,5);
trial_label=ActionTrain(:,1);
trial_length=min(end_feedback-start_feedback);

freq=zeros(size(Final_features,1),1);

for d=1:size(Final_features,1)
 freq(d)=find(Frequencies==Final_features(d,1));
 end

selected_features=[freq,Final_features(:,2)];

%%
all_features=[];%zeros(length(start_feedback),trial_length,size(selected_features,1)); % trial x window x features
all_label=[];%zeros(length(start_feedback),trial_length,2);

for i=1:length(start_feedback)
    
    Window_single_feat=[];
    
    %j'extrais tout les features pour une cue
    for a=1:size(selected_features,1)
    Window_single_feat=[Window_single_feat,selected_data(start_feedback(i):end_feedback(i),selected_features(a,1),selected_features(a,2))];
    end
    
    
    %je fais une matrice 3D [trial x time x features] donc chaque trial
    %doit avoir la meme duree
    trans=Window_single_feat(1:trial_length,:);
    
    all_features= [all_features;trans];
    
    %je fais une matrice 3D [trial x time x cue_label|trial_label]
    all_label=[all_label;ones(trial_length,1).*ActionTrain(i,1)];
 end
    
%%
    classifier = fitcdiscr(all_features, all_label, 'discrimtype', 'linear'); %train an LDA classifier

    save('Classifier\Final classifier','classifier');
    Param=[OptimalAlpha(1),OptimalLim(1)]
    Features_index=selected_features;
    
    save('Classifier\Final classifier Param','Param');

    save('Classifier\Final classifier Features','Features_index');


    
function [class_error]=error(real_label, predicted_label)
    false=nnz(real_label-predicted_label);
    
    class_error=(false/length(real_label));
end


