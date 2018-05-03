function GetPSD_TimePlotTry(Epoch_both_feet, Epoch_both_hands, Baseline_both_feet, Baseline_both_hands)

tot_num_elec = size(Epoch_both_feet,4); 
figure()
for n_elec = 1:tot_num_elec
    subplot(4,4,n_elec)
    
    % average the psd over frequencies and over trials
    PSD_both_feet = mean(mean(Epoch_both_feet(:,:,:,n_elec),3),2)
    PSD_both_hands = mean(mean(Epoch_both_hands(:,:,:,n_elec),3),2); 
    
    PSD_baseline_both_feet = mean(mean(Baseline_both_feet(:,:,:,n_elec),3),2);
    PSD_baseline_both_hands = mean(mean(Baseline_both_hands(:,:,:,n_elec),3),2);
    
    % plot
    plot(log10([PSD_baseline_both_feet; PSD_both_feet])); 
    hold on;
    plot(log10([PSD_baseline_both_hands; PSD_both_hands]));
    hold on;
    fixation_duration = length(PSD_baseline_both_hands)
    line([fixation_duration fixation_duration], get(gca, 'ylim'),'Color','green','LineStyle','--')
    
    title(['Electrode ', num2str(n_elec)])  
end
legend('both feet','both hands','end of fixation');


end

