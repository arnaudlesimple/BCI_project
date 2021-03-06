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

%%

%Je recupere les temps de chaque trial ainsi que son label
start_feedback=Action(:,4);
end_feedback=Action(:,5);
trial_label=Action(:,1);
trial_length=min(end_feedback-start_feedback);

% en fonction de la discriminancy je choisis les features [features x freq electrode]
 features=[24 9]; %[frequ x channel]
 freq=zeros(size(features,1),1);
 
 for d=1:size(features,1)
 freq(d)=find(Frequencies==features(d,1));
 end
% selected_features=[selected_features(:,1)./2-1,selected_features(:,2)]; %[frequ_INDEX x channel

selected_features=[freq,features(:,2)];
% for a=1:size(features,1)
%     freq(a)=find(Frequencies==features(a,1));
% end
%     selected_features=[freq,features(:,2)];

%%
window_feat=zeros(length(start_feedback),trial_length,size(selected_features,1)); % trial x window x features
window_label=zeros(length(start_feedback),trial_length,2);

for i=1:length(start_feedback)
    
    Window_single_feat=[];
    
    %j'extrais tout les features pour une cue
    for a=1:size(selected_features,1)
    Window_single_feat=[Window_single_feat,selected_data(start_feedback(i):end_feedback(i),selected_features(a,1),selected_features(a,2))];
    end
    
    
    %je fais une matrice 3D [trial x time x features] donc chaque trial
    %doit avoir la meme duree
    
    window_feat(i,:,:)=Window_single_feat(1:trial_length,:);
    
    %je fais une matrice 3D [trial x time x cue_label|trial_label]
    window_label(i,:,:)=ones(trial_length,2).*[Action(i,1),i];
 end

%%
%partition sur les trials afin d'avoir toutes les windows dans le meme set
K=10;
partition = cvpartition(length(start_feedback), 'KFold', K);


%%
single_sample_accuracy=zeros(K,1);
trial_accuracy=zeros(K,1);

for i=1:K
    %je veux avoir les windows des trials ensemble dans les sets
    
    train_trials=find(partition.training(i));
    test_trials=find(partition.test(i)); %donne index 1 de la matrice 3D
    
    training_set=[];
    testing_set=[];
   
    train_label=[];
    test_label=[];
    
    % pour utiliser fitcdiscr il faut une 2D donc je mets � la suite
    % label =[time x cue_label|trial_label]
    
    for n=1:partition.TrainSize(i)
        training_set=[training_set,squeeze(window_feat(train_trials(n),:,:))];
        train_label=[train_label;squeeze(window_label(train_trials(n),:,:))];
    end
    
    for n=1:partition.TestSize(i)
        testing_set=[testing_set;squeeze(window_feat(test_trials(n),:,:))];
        test_label=[test_label;squeeze(window_label(test_trials(n),:,:))];
    end
    
    classifier = fitcdiscr(training_set', train_label(:,1), 'discrimtype', 'linear'); %train an LDA classifier

    [predicted_label,score,cost] = predict(classifier, testing_set); %score [771 773]

    single_sample_accuracy(i) =  accuracy( test_label(:,1), predicted_label);

    %trial accuracy
    alpha=0.9;
    limit771=0.8;
    limit773=0.2;
    
    decision771=zeros(trial_length,length(test_trials)); %proba     [time x trial]

    final_classification=zeros(length(test_trials),1); %label apr�s estimation de chaque trial
    TestTrialLabel=window_label(test_trials,1,1);
    
        for m=1:length(test_trials)

            trial_number=test_trials(m); %donne indice 1 de la matrice 3D            
            ProbaTrial=score(test_label(:,2)==trial_number,1);
            
            for l=2:length(ProbaTrial)

                decision771(l,m)=ProbaTrial(l)*alpha+decision771(l-1,m)*(1-alpha);

                if decision771(l,m) > limit771

                    final_classification(m)=771;
                else if decision771(l,m) < limit773
                      final_classification(m)=773;  
                    end
                end
            end
        end
        
     trial_accuracy(i,1)=accuracy(TestTrialLabel,final_classification)
     
    
end

     avg_s_s_accuracy=mean(single_sample_accuracy)
     avg_trial_accuracy=mean(trial_accuracy)



%%

plot(decision771(:,1))
hold on;
plot(decision771(:,2))
hold on;
plot(decision771(:,3))
hold on;
plot([1:trial_length],ones(trial_length)*limit771)
hold on;
plot([1:trial_length],ones(trial_length)*limit773)

%%

function [class_accuracy]=accuracy(real_label, predicted_label)
    false=nnz(real_label-predicted_label);
    
    class_accuracy=1-(false/length(real_label));
end
    