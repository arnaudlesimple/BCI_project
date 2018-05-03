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

mu_band = [3:6];
beta_band = [7:18];
mu_beta_band = [3:18];
all_band=1:23;

band = {mu_band, beta_band};
band_selected = all_band;


%%

discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies);

%%




start_feedback=Action(:,4);
end_feedback=Action(:,5);

trial_label=Action(:,1);

selected_features= [24 9; 12 2; 12 3; 12 5; 12 7;12 8;12 11]; %[frequ x channel]

selected_features=[selected_features(:,1)./2-1,selected_features(:,2)]; %[frequ_INDEX x channel


%%
window_feat=[]; % window x features
window_label=[];

for i=1:length(start_feedback)
    
    Window_single_feat=[];
    
    for a=1:length(selected_features)
    Window_single_feat=[Window_single_feat,selected_data(start_feedback(i):end_feedback(i),selected_features(a,2),selected_features(a,1))];
    end
    
    window_feat=[window_feat;Window_single_feat];
    
   temp=ones(end_feedback(i)-start_feedback(i)+1,2).*[Action(i,1),i]

    
    window_label=[window_label; temp];
 end

%%

partition_N = cvpartition(length(start_feedback), 'KFold', 10);


erreur=zeros(10,1);

for i=1:10
%on recherche les index des training/test sample
training_set_index = find(partition_N.training(i));
testing_set_index = find(partition_N.test(i));

features_training=window_feat(window_label(:,2)==;
features_testing=window_feat(testing_set_index,:);

classifier = fitcdiscr(features_training, window_label(training_set_index), 'discrimtype', 'linear'); %train an LDA classifier
test_label = predict(classifier, features_testing); 

single_sample_accuracy(i) =  accuracy( window_label(testing_set_index), test_window_label);




end

avg_s_s_accuracy=mean(single_sample_accuracy);





function [class_accuracy]=accuracy(real_label, predicted_label)
    false=nnz(real_label-predicted_label);
    
    class_accuracy=1-(false/length(real_label));
end