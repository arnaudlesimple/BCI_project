clear
close all


Action = load('SPD/Event Window.mat');
Action = Action.Event_window;
Event=load('SPD/WindowLabel.mat');%WindowLabel
Event=Event.labelAction;

%% Attention alpha doit etre penser en terme de practicité car sinon on a pas le temps
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
    
            [mean_window(Classifier,i,j),std_window,mean_trial(Classifier,i,j),std_trial]=GetClassifierAccuracy(selected_data,Action,features,alpha(j),1-limit(i),limit(i),Frequencies,Classifier) ;
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
[a1,b1]=min(mean_trial(1,:,:));
[a2,b2]=min(min(mean_trial(1,:,:)));
a=a1(b2);
b=b1(b2); %optimal_limit=b et b2=optimal alpha
%surf([limit(b)-0.1:0.1:limit(b)+0.1],[alpha(b2)-0.1:0.1:alpha(b2)+0.1],ones(3,3)*A(b,b2));
surf([limit(b)-0.025,limit(b),limit(b)+0.025],[alpha(b2)-0.05,alpha(b2),alpha(b2)+0.05],ones(3,3)*A(b,b2));

hold on;
surf(limit,alpha,A);
optimal=['minimal error is ' num2str(a) ' and optimal limit =' num2str(limit(b)) 'with alpha= ' num2str(alpha(b2))];
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

%% ffs feature



%% FFS

alpha=0.1;
limit771=0.7;
limit773=0.3;

Feat_selected=zeros(4,size(features,1));
mean_trial=zeros(4,1);

for Classifier=1:4
     [Feat_selected(Classifier,:),historyCrit,historyIn] = FFS(selected_data,Action,features,OptimalAlpha(Classifier),1-OptimalLim(Classifier),OptimalLim(Classifier),Frequencies,Classifier);
     mean_trial(Classifier)=historyCrit(end);
end

%% d
Final= [mean_trial,Feat_selected]
Final_features_index=find(Feat_selected(a,:));
Final_features=features(Final_features_index,:)

%% build final classifier
    a=1;
     
    
    %%
    %Je recupere les temps de chaque trial ainsi que son label
start_feedback=Action(:,4);
end_feedback=Action(:,5);
trial_label=Action(:,1);
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
    all_label=[all_label;ones(trial_length,1).*Action(i,1)];
 end
    
%%
    classifier = fitcdiscr(all_features, all_label, 'discrimtype', 'linear'); %train an LDA classifier

    save('Classifier\Final classifier','classifier');
    Param=[OptimalAlpha(1),OptimalLim(1)]
    Features_index=selected_features;
    
    save('Classifier\Final classifier Param','Param');

    save('Classifier\Final classifier Features','Features_index');




