function [avg_s_s_accuracy,STD_s_s_accuracy,avg_trial_accuracy,STD_trial_accuracy,mean_Confusion,AUC] = GetClassifierAccuracy(selected_data,Action,features,alpha,limit771,limit773,Frequencies_set,Classifier)


switch Classifier
    case 1
        Classifier='linear'
    case 2
        Classifier='diaglinear'
    case 3
        Classifier='quadratic'
    case 4
        Classifier='diagquadratic'
end

%Je recupere les temps de chaque trial ainsi que son label
start_feedback=Action(:,4);
end_feedback=Action(:,5);
trial_label=Action(:,1);
trial_length=min(end_feedback-start_feedback);


freq=zeros(size(features,1),1);

for d=1:size(features,1)
 freq(d)=find(Frequencies_set==features(d,1));
 end

selected_features=[freq,features(:,2)];

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
Confusion_M=zeros(K,3,3);

AUCwindow=[];

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
    
    if size(features,1) ~=1
        for n=1:partition.TrainSize(i)
            training_set=[training_set;squeeze(window_feat(train_trials(n),:,:))];
            train_label=[train_label;squeeze(window_label(train_trials(n),:,:))];
        end

        for n=1:partition.TestSize(i)
            testing_set=[testing_set;squeeze(window_feat(test_trials(n),:,:))];
            test_label=[test_label;squeeze(window_label(test_trials(n),:,:))];
        end
    end
    
    %car window_feat a desormais 2D
    if size(features,1) ==1
        
        for n=1:partition.TrainSize(i)
            test=squeeze(window_feat(train_trials(n),:));
            test2=window_feat(train_trials(n),:);
            training_set=[training_set;window_feat(train_trials(n),:)'];
            train_label=[train_label;squeeze(window_label(train_trials(n),:,:))];
        end

        for n=1:partition.TestSize(i)
            testing_set=[testing_set;window_feat(test_trials(n),:)'];
            test_label=[test_label;squeeze(window_label(test_trials(n),:,:))];
        end  
    end
    
    classifier = fitcdiscr(training_set, train_label(:,1), 'discrimtype', Classifier); %train an LDA classifier

    [Predicted_label,score,cost] = predict(classifier, testing_set); %score [771 773]
    
    
    if sum(test_label(:,1)==771) ~=0 && sum(test_label(:,1)==773) ~=0
        
    [Xt,Yt,T,AUCt]=perfcurve(test_label(:,1),score(:,1),'771');

     AUCwindow=[AUCwindow,AUCt];
    end

    single_sample_accuracy(i) =  error( test_label(:,1), Predicted_label);

    %trial accuracy
    
    %!!!!! commence a 0.5
    %decision771=zeros(trial_length,length(test_trials)); %proba     [time x trial]
    decision771=ones(trial_length,length(test_trials))*0.5; %proba     [time x trial]


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
        
     trial_accuracy(i,1)=error(TestTrialLabel,final_classification);

     Confusion_M(i,:,:) = confusionmat(TestTrialLabel,final_classification,'Order',[0 771 773]);
     
    
end

     avg_s_s_accuracy=mean(single_sample_accuracy);
     avg_trial_accuracy=mean(trial_accuracy);
     
     STD_s_s_accuracy=std(single_sample_accuracy);
     STD_trial_accuracy=std(trial_accuracy);
     
     mean_Confusion=squeeze(mean(Confusion_M,1));
     

    AUC=mean(AUCwindow);




function [class_error]=error(real_label, predicted_label)
    false=nnz(real_label-predicted_label);
    
    class_error=(false/length(real_label));
end

end

