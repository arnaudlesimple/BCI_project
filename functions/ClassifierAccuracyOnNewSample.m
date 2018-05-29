function [avg_s_s_error,STD_s_s_error,avg_trial_error,STD_trial_error,FINAL,Confusion_M,Xwindow,Ywindow,AUCwindow] = ClassifierAccuracyOnNewSample(DataTest,DataTrain,ActionTest,ActionTrain,features,alpha,limit771,limit773,Frequencies_set,Classifier)


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
start_feedback_Train=ActionTrain(:,4);
end_feedback_Train=ActionTrain(:,5);
trial_label_Train=ActionTrain(:,1);
trial_length_Train=min(end_feedback_Train-start_feedback_Train);

start_feedback_Test=ActionTest(:,4);
end_feedback_Test=ActionTest(:,5);
trial_label_Test=ActionTest(:,1);
trial_length_Test=min(end_feedback_Test-start_feedback_Test);

freq=zeros(size(features,1),1);

for d=1:size(features,1)
 freq(d)=find(Frequencies_set==features(d,1));
 end
% selected_features=[selected_features(:,1)./2-1,selected_features(:,2)]; %[frequ_INDEX x channel

selected_features=[freq,features(:,2)];

%%
window_feat_Train=zeros(length(start_feedback_Train),trial_length_Train,size(selected_features,1)); % trial x window x features
window_label_Train=zeros(length(start_feedback_Train),trial_length_Train,2);

window_feat_Test=zeros(length(start_feedback_Test),trial_length_Test,size(selected_features,1)); % trial x window x features
window_label_Test=zeros(length(start_feedback_Test),trial_length_Test,2);

for i=1:length(start_feedback_Train)
    
    Window_single_feat_Train=[];
    
    %j'extrais tout les features pour une cue
    for a=1:size(selected_features,1)
    Window_single_feat_Train=[Window_single_feat_Train,DataTrain(start_feedback_Train(i):end_feedback_Train(i),selected_features(a,1),selected_features(a,2))];
    end
    
    
    %je fais une matrice 3D [trial x time x features] donc chaque trial
    %doit avoir la meme duree
    
    window_feat_Train(i,:,:)=Window_single_feat_Train(1:trial_length_Train,:);
    
    %je fais une matrice 3D [trial x time x cue_label|trial_label]
    window_label_Train(i,:,:)=ones(trial_length_Train,2).*[ActionTrain(i,1),i];
end

 for i=1:length(start_feedback_Test)
    
    Window_single_feat_Test=[];
    
    %j'extrais tout les features pour une cue
    for a=1:size(selected_features,1)
    Window_single_feat_Test=[Window_single_feat_Test,DataTest(start_feedback_Test(i):end_feedback_Test(i),selected_features(a,1),selected_features(a,2))];
    end
    
    
    %je fais une matrice 3D [trial x time x features] donc chaque trial
    %doit avoir la meme duree
    
    window_feat_Test(i,:,:)=Window_single_feat_Test(1:trial_length_Test,:);
    
    %je fais une matrice 3D [trial x time x cue_label|trial_label]
    window_label_Test(i,:,:)=ones(trial_length_Test,2).*[ActionTest(i,1),i];
 end
%%
%partition sur les trials afin d'avoir toutes les windows dans le meme set



%%
%     train_trials=1:60;
%     test_trials=
    
    
single_sample_error=zeros(1,1);
trial_error=zeros(1,1);
    
    training_set=[];
    testing_set=[];
   
    train_label=[];
    test_label=[];
    
    % pour utiliser fitcdiscr il faut une 2D donc je mets à la suite
    % label =[time x cue_label|trial_label]
    
    if size(features,1) ~=1
        for n=1:length(start_feedback_Train)
            training_set=[training_set;squeeze(window_feat_Train(n,:,:))];
            train_label=[train_label;squeeze(window_label_Train(n,:,:))];
        end

        for n=1:length(start_feedback_Test)
            testing_set=[testing_set;squeeze(window_feat_Test(n,:,:))];
            test_label=[test_label;squeeze(window_label_Test(n,:,:))];
        end
    end
    
    %car window_feat a desormais 2D
    if size(features,1) ==1
        
        for n=1:length(start_feedback_Train)
            test=squeeze(window_feat_Train(n,:));
            test2=window_feat_Train(n,:);
            training_set=[training_set;window_feat_Train(n,:)'];
            train_label=[train_label;squeeze(window_label_Train(n,:,:))];
        end

        for n=1:length(start_feedback_Test)
            testing_set=[testing_set;window_feat_Train(n,:)'];
            test_label=[test_label;squeeze(window_label_Train(n,:,:))];
        end  
    end
    
    classifier = fitcdiscr(training_set, train_label(:,1), 'discrimtype', Classifier); %train an LDA classifier
    
    FINAL=fitcdiscr([training_set;testing_set], [train_label(:,1);test_label(:,1)], 'discrimtype', Classifier);
    

    [Predicted_label,score,cost] = predict(classifier, testing_set); %score [771 773]
    
    [Xwindow,Ywindow,T,AUCwindow]=perfcurve(test_label(:,1),score(:,1),'771');
    

    
    
    single_sample_error =  error( test_label(:,1), Predicted_label);

    %trial accuracy
    trial_length=size(window_feat_Test,2);
    decision771=ones(trial_length,size(window_feat_Test,1))*0.5; %proba     [time x trial]

    final_classification=zeros(size(window_feat_Test,1),1); %label après estimation de chaque trial
    TestTrialLabel=trial_label_Test;
    
        for m=1:size(window_feat_Test,1)
            
            ProbaTrial=score(test_label(:,2)==m,1);
            
            %trial_number=test_trials(m); %donne indice 1 de la matrice 3D            
            
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
        
     trial_error=error(TestTrialLabel,final_classification);
     
    Confusion_M = confusionmat(TestTrialLabel,final_classification,'Order',[0 771 773]);


     avg_s_s_error=mean(single_sample_error);
     avg_trial_error=mean(trial_error);
     
     STD_s_s_error=std(single_sample_error);
     STD_trial_error=std(trial_error);
    




function [class_error]=error(real_label, predicted_label)
    false=nnz(real_label-predicted_label);
    
    class_error=(false/length(real_label));
end

end

