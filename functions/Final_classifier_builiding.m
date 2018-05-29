clear
close all
clc

%we train on offline data, measured the same day

name='Elise';

addpath(['C:\Users\Mathieu\Documents\EPFL\Master\MA2\BCI\PSD final\' name ])
% offline: 1trial = 30 cues
load('offline\Event Window.mat');
OfflineAction=Event_window;
load('offline\WindowLabelDay.mat');
OfflineDayLabel=WindowLabelDay;
load('offline\WindowLabel.mat');
OfflineActionLabel=labelAction;
load('offline\WindowLabelRun.mat');
OfflineRunLabel=WindowLabelRun;
load('offline\SPD with SmallLaplacian Spatial filtre.mat');
OfflineData=psdt;

% online: 1 trial =20 cues
load(['C:\Users\Mathieu\Documents\EPFL\Master\MA2\BCI\PSD final\' name '\online\Event Window.mat']);
OnlineAction=Event_window;
load(['C:\Users\Mathieu\Documents\EPFL\Master\MA2\BCI\PSD final\' name '\online\WindowLabelDay.mat']);
OnlineDayLabel=WindowLabelDay;
load(['C:\Users\Mathieu\Documents\EPFL\Master\MA2\BCI\PSD final\' name '\online\WindowLabel.mat']);
OnlineActionLabel=labelAction;
load(['C:\Users\Mathieu\Documents\EPFL\Master\MA2\BCI\PSD final\' name '\online\WindowLabelRun.mat']);
OnlineRunLabel=WindowLabelRun;
load(['C:\Users\Mathieu\Documents\EPFL\Master\MA2\BCI\PSD final\' name '\online\SPD with SmallLaplacian Spatial filtre.mat']);
OnlineData=psdt;


window_frequency = 16;
load('SPD/Anonymous/Frequences.mat');
frequencies = load('SPD/Anonymous/Frequences.mat');


mu_band = 3:6;
beta_band = 7:18;
mu_beta_band = 3:18;
all_band=1:23;

band = {mu_band, beta_band};
band_selected = all_band;

Classifier=[1,2,3,4]; %["linear","diaglinear","quadratic","diagquadratic"];
%% si tu veux mettre du 
TotalTrial=size(OfflineRunLabel,1)+size(OnlineRunLabel,1);
NumberTrain=round(0.7*TotalTrial)
NumberTest=TotalTrial-NumberTrain

% To modify yourself
ActionTrain=OfflineAction;
ActionTest=OnlineAction;


%% Observation of Train discriminancy map
TrialNum=size(ActionTrain,1)/30;
for a= 1:TrialNum

discrimancy = GetDiscrimancyMap(OfflineData, band_selected, window_frequency, frequencies,'inutile',ActionTrain((a-1)*30+1:a*30,:));
end
%% Discr for all
discrimancy = GetDiscrimancyMap(OfflineData, band_selected, window_frequency, frequencies,'inutile',ActionTrain);

%% Mathieu
features= [22 13; 24 15; 32 13; 30 13; 8 13; 22 6; 30 7; 6 7]; %[frequ x channel]
%% ELise
features=  [12 7; 10 7;8 11; 12 11; 10 11; 8 7];
%% pour que ce soit dans le meme ordre que après
features=sortrows(features,2);

alpha=0.06;
limit771=0.8;
limit773=0.2;
mean_trial=zeros(4,1);

Feat_selected=zeros(4,length(features));
%% test 4 classifier
%% FFS features

for Classifier=1:4
     [Feat_selected(Classifier,:),historyCrit,historyIn] = FFS(OfflineData,ActionTrain,features,alpha,limit771,limit773,Frequencies,Classifier);
     mean_trial(Classifier)=historyCrit(end);
end

%% d
Final= [mean_trial,Feat_selected]


%% print Final

ClassifierName = {'Linear';'Diaglinear';'Quadratic';'Diagquadratic'};
Validation_error=round(Final(:,1),3);
F1=Feat_selected(:,1);
F2=Feat_selected(:,2);
F3=Feat_selected(:,3);
F4=Feat_selected(:,4);
F5=Feat_selected(:,5);
F6=Feat_selected(:,6);
N_feature_chosen=sum(Feat_selected,2);


T = table(Validation_error,F1,F2,F3,F4,F5,F6,N_feature_chosen,'RowNames',ClassifierName)%
T2 = table(Validation_error,N_feature_chosen,'RowNames',ClassifierName)%


%%

uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
%%
uitable('Data',T2{:,:},'ColumnName',T2.Properties.VariableNames,...
    'RowName',T2.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

%%
FFS_features={features(find(Feat_selected(1,:)),:),features(find(Feat_selected(2,:)),:),features(find(Feat_selected(3,:)),:),features(find(Feat_selected(4,:)),:)};


%% different alpha

alpha=0.01:0.01:0.07;
mean_window=zeros(length(alpha),4);
mean_trial=zeros(length(alpha),4);
mean_confusion=zeros(length(alpha),3,3);
AUC_alpha=zeros(length(alpha),4);



for Classifier=1:4


    for i=1:length(alpha)

        [mean_window(i,Classifier),std_window,mean_trial(i,Classifier),std_trial,mean_confusion(i,:,:),AUC_alpha(i)]=GetClassifierAccuracy(OfflineData,ActionTrain,FFS_features{Classifier},alpha(i),limit771,limit773,Frequencies,Classifier);
    end
            [a,b]=min(mean_trial(:,Classifier));
            
            best_Confusion(Classifier,:,:)=mean_confusion(b,:,:);
            AUC_class(Classifier)=AUC_alpha(b);

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
ylabel('Validation error [%]')

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
ylabel('Validation error [%]')

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
ylabel('Validation error [%]')

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
ylabel('Validation error [%]')

%% confusion matrix and AUC
%for Class=1:4



Conf=squeeze(best_Confusion(1,:,:));
Cue={'Over Time';'Hand cue';'Feet cue'}
Over=round(Conf(:,1));
Hand=round(Conf(:,2));
Both_Feet=round(Conf(:,3));
T = table(Over,Hand,Both_Feet,'RowNames',Cue)%

hf = figure;
ha = subplot(2,2,1);
pos = get(ha,'Position');
un = get(ha,'Units');
axis off
A=uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',pos);
A.FontSize=12
titre=['Linear classifier, AUC=' num2str(round(AUC_class(1),3))];
title(titre)


Conf=squeeze(best_Confusion(2,:,:));
Cue={'Over Time';'Hand cue';'Feet cue'}
Over=round(Conf(:,1));
Hand=round(Conf(:,2));
Both_Feet=round(Conf(:,3));
T = table(Over,Hand,Both_Feet,'RowNames',Cue)%
%hf = figure;
ha = subplot(2,2,2);
pos = get(ha,'Position');
un = get(ha,'Units');
axis off
A=uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',pos);
A.FontSize=12
titre=['Diaglinear classifier, AUC=' num2str(round(AUC_class(2),3))];
title(titre)



Conf=squeeze(best_Confusion(3,:,:));
Cue={'Over Time';'Hand cue';'Feet cue'}
Over=round(Conf(:,1));
Hand=round(Conf(:,2));
Both_Feet=round(Conf(:,3));
T = table(Over,Hand,Both_Feet,'RowNames',Cue)%
%hf = figure;
ha = subplot(2,2,3);
pos = get(ha,'Position');
un = get(ha,'Units');
axis off
A=uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',pos);
A.FontSize=12
titre=['Quadratic classifier, AUC=' num2str(round(AUC_class(3),3))];
title(titre)



Conf=squeeze(best_Confusion(4,:,:));
Cue={'Over Time';'Hand cue';'Feet cue'}
Over=round(Conf(:,1));
Hand=round(Conf(:,2));
Both_Feet=round(Conf(:,3));
T = table(Over,Hand,Both_Feet,'RowNames',Cue)%
%hf = figure;
ha = subplot(2,2,4);
pos = get(ha,'Position');
un = get(ha,'Units');
axis off
A=uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',pos);
A.FontSize=12
titre=['Diagquadratic classifier, AUC=' num2str(round(AUC_class(4),3))];
title(titre)




%% test on new data
Param={OptimalAlpha;FFS_features}

Classifier=3;


[avg_s_s_error,STD_s_s_error,avg_trial_error,STD_trial_error,FINALClassifier,Confusion_M,Xwindow,Ywindow,AUCwindow]=ClassifierAccuracyOnNewSample(OnlineData,OfflineData,ActionTest,ActionTrain,FFS_features{Classifier},OptimalAlpha(Classifier),limit771,limit773,Frequencies,Classifier)

%% image final test

    plot(Xwindow,Ywindow);   

    ylabel('True positive rate');
    xlabel('False negative rate');
    
    leg=['ROC curve of final classifier, AUC=' num2str(AUCwindow) '.']
    legend(leg);
    
    
     %title('Final classifier ROC curve (window prediction)');

    
 %% final conf
 
 Conf=Confusion_M;
Cue={'Hand cue';'Feet cue'}
Unconclusive=round(Conf(2:3,1));
Both_Hand=round(Conf(2:3,2));
Both_Feet=round(Conf(2:3,3));
T = table(Unconclusive,Both_Hand,Both_Feet,'RowNames',Cue)%

A=uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

A.FontSize=13;

%%  final test error


Cue={'Final classifier test error (Quadratic)'}
Window_test_error=round(avg_s_s_error,2);
Trial_test_error=round(avg_trial_error,2);

T = table(Window_test_error,Trial_test_error,'RowNames',Cue)%

B=uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
B.FontSize=13;

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


