function getFrequencyTimePlot(psd, band, window_frequency, frequencies)
    
    % Extract epoch 
    [Epoch_both_feet, Epoch_both_hands, Baseline_both_feet, Baseline_both_hands, trial_length_feet, trial_length_hand] = Epoching(psd, band, files);
    number_electrode = 16;

    % Both feet: for each electrode
    figure()
    for n_electrode = 1:number_electrode
        PSD_both_feet = mean(Epoch_both_feet(:,:,:,n_electrode),3);
        BaseLine = mean(mean(Baseline_both_feet(:,:,:,n_electrode),3));
        ERD_ERS = 10.*(PSD_both_feet - BaseLine) ./ BaseLine;
        subplot(4,4,n_electrode)
        imagesc([0:trial_length_feet/window_frequency],frequencies.Frequencies(band), transpose(ERD_ERS))
        c = colorbar; 
        %c.Limits = [-0.5, 0.5]
        maxi(n_electrode) = max(max(ERD_ERS));
        mini(n_electrode) = min(min(ERD_ERS));
        xlabel('time [s]')
        ylabel('frequency [Hz]')
        title(['electrode' mat2str(n_electrode)])
        
        % find limits of colorbar to have the same legend for each electrode
        set(gca, 'XTick', [0:trial_length_feet/window_frequency])
        set(gca, 'YTick', frequencies.Frequencies(band), 'fontsize', 5)
    end
    [min(mini),max(maxi)]

    % Both hand: for each electrode
    figure()
    for n_electrode = 1:number_electrode
        PSD_both_hand = mean(Epoch_both_hands(:,:,:,n_electrode),3);
        BaseLine = mean(mean(Baseline_both_hands(:,:,:,n_electrode),3));
        ERD_ERS = 10.*(PSD_both_hand - BaseLine) ./ BaseLine;

        subplot(4,4,n_electrode)
        imagesc([0:trial_length_hand/window_frequency],frequencies.Frequencies(band), transpose(ERD_ERS))
        c = colorbar;
        %c.Limits = [-0.5,0.5];
        xlabel('time [s]')
        ylabel('frequency [Hz]')
        set(gca, 'XTick', [0:trial_length_hand/window_frequency])
        set(gca, 'YTick', frequencies.Frequencies(band), 'fontsize', 5)
        title(['electrode' mat2str(n_electrode)])
        
        % find limits of colorbar to have the same legend for each electrode
        maxi(n_electrode) = max(max(ERD_ERS));
        mini(n_electrode) = min(min(ERD_ERS));
    end
    [min(mini),max(maxi)]

end