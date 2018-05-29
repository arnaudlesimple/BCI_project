function discrimancy = GetDiscrimancyMap(psd, band, window_frequency, frequencies, files, Event_window)

    % Extract epoch 
    [Epoch_both_feet, Epoch_both_hands, Baseline_both_feet, Baseline_both_hands, trial_length_feet, trial_length_hand] = Epoching(psd, band, Event_window);
    number_electrode = 16;
        
    for n_electrode = 1:number_electrode
        % Select PSD for one channel with different windows and frequencies
        PSD_both_feet = mean(Epoch_both_feet(:,:,:,n_electrode),3); % mean along all trials
        BaseLine_both_feet = mean(mean(Baseline_both_feet(:,:,:,n_electrode),3));
        ERD_ERS_both_feet = 10.*(PSD_both_feet - BaseLine_both_feet)./ BaseLine_both_feet;
        
        PSD_both_hands = mean(Epoch_both_hands(:,:,:,n_electrode),3); % mean along all trials
        BaseLine_both_hands = mean(mean(Baseline_both_hands(:,:,:,n_electrode),3));
        ERD_ERS_both_hands = 10.*(PSD_both_hands - BaseLine_both_hands) ./ BaseLine_both_hands;
        
        for f = 1:size(PSD_both_feet,2)
            %plot(ERD_ERS_both_feet(:,f))
            %plot(ERD_ERS_both_hands(:,f))                       
            discrimancy(f, n_electrode) = abs(mean(ERD_ERS_both_feet(:,f)) - mean(ERD_ERS_both_hands(:,f))) / sqrt(std(ERD_ERS_both_feet(:,f))^2 + std(ERD_ERS_both_hands(:,f))^2);
        end
        
%         figure()
%         plot(mean(ERD_ERS_both_feet))
%         hold on
%         plot(mean(ERD_ERS_both_hands))
%         xlabel('time')
%         ylabel('PSD')
    end
    
    figure()
    imagesc([1:number_electrode],frequencies.Frequencies(band), discrimancy)
    colorbar
    xlabel('channel')
    ylabel('frequency')
    title('Discriminancy map');
    set(gca, 'YTick', frequencies.Frequencies(band), 'fontsize', 5)
    %set(gca, 'XTick', frequencies.Frequencies(band), 'fontsize', 5)
end

