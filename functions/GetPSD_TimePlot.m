function GetPSD_TimePlot(psdt, labelAction, Frequencies, band_freq, window_label, files, Event_window)

%% parameters
% load(['SPD/' files '/Event Window.mat']);
% % Take only windows that correspond to the selected run
% Event_window = Event_window(find((Event_window(:,2) > run_ind(1)) & (Event_window(:,5) < run_ind(end))),:);

%table with all both_feet trials
index_trial_feet = find(Event_window==771);
num_trial_feet = length(index_trial_feet);
trial_feet = zeros(num_trial_feet, size(Event_window,2)-1);
trial_feet_fixation = zeros(num_trial_feet, 1);

for i = 1:num_trial_feet
    trial_feet(i,:)= Event_window(index_trial_feet(i), 2:size(Event_window,2));
    trial_feet_fixation(i) = Event_window(index_trial_feet(i), 3)-Event_window(index_trial_feet(i), 2);
end

feet_mean_fixation_duration = mean(trial_feet_fixation);
%table with all both_hands trials
index_trial_hands = find(Event_window==773);
num_trial_hands = length(index_trial_hands);
trial_hands = zeros(num_trial_hands, size(Event_window,2)-1);
for i = 1:num_trial_hands
    trial_hands(i,:)= Event_window(index_trial_hands(i), 2:size(Event_window,2));
    trial_hands_fixation(i) = Event_window(index_trial_hands(i), 3)-Event_window(index_trial_hands(i), 2);
end
hands_mean_fixation_duration = mean(trial_hands_fixation);


%minimal duration among all trials for each event
trial_length_feet = min(trial_feet(:,size(trial_feet,2))-trial_feet(:,1)); 
trial_length_hands = min(trial_hands(:,size(trial_hands,2))-trial_hands(:,1)); 
%minimal duration among all trials for the fixation
fixation_duration=mean([hands_mean_fixation_duration,feet_mean_fixation_duration]);
%% Get PSD values for each event

tot_num_elec = size(psdt,3);
figure()
for n_electrode=1:tot_num_elec

    subplot(4,4,n_electrode)

    PSD_both_feet = zeros(trial_length_feet,length(band_freq),num_trial_feet); %PDS vs time for each band and for each trial
    PSD_both_hands = zeros(trial_length_hands,length(band_freq),num_trial_hands);

    % Epoch_both_hands = ;
     Epoch_both_feet = zeros(trial_length_feet,length(band_freq),num_trial_feet);

    for trial_number = 1:num_trial_feet
        %feet_trial_samples = [trial_feet(trial_number,1):trial_feet(trial_number,size(trial_feet,2))];
        % Extract mean PSD over time from both feet
        Epoch_both_feet(:,:,trial_number) = psdt(trial_feet(trial_number,1):(trial_feet(trial_number,1)+trial_length_feet-1), band_freq ,n_electrode);

    end

    for trial_number = 1:num_trial_hands
        Epoch_both_hands(:,:,trial_number) = psdt(trial_hands(trial_number,1):(trial_hands(trial_number,1)+trial_length_hands-1), band_freq ,n_electrode);
    end

    %% average the psd over frequencies and over trials

    PSD_both_feet = mean(mean(Epoch_both_feet(:,:,:),3),2)
    PSD_both_hands = mean(mean(Epoch_both_hands(:,:,:),3),2);

    %% plot
    plot(log10(PSD_both_feet)); hold on;
    plot(log10(PSD_both_hands));
    hold on;
    line([fixation_duration fixation_duration], get(gca, 'ylim'),'Color','green','LineStyle','--')
    title(['Electrode ', num2str(n_electrode)])
end
legend('both feet','both hands','end of fixation');

end