function GetFeaturePlot(psd, band, frequencies, window_label, psd_num, band_num)
    
    psd_type = {'small laplacian', 'large laplacian', 'CAR', 'no spatial filter'}
    band_type = {'mu band', 'beta band'};

    % Extract mean PSD over time from both feet
    PSD_both_feet = mean(psd.psdt(window_label.labelAction==771, band, :),1);
    % Extract mean PSD over time from both hand
    PSD_both_hand = mean(psd.psdt(window_label.labelAction==773, band, :),1);
    % reshape PSD matrices to go from 3D to 2D in order to use imagesc
    S = size(PSD_both_feet);
    PSD_both_feet = reshape(PSD_both_feet,[S(1)*S(2),S(3)]); 
    PSD_both_hand = reshape(PSD_both_hand,[S(1)*S(2),S(3)]); 
    
    % Compute the ratio of power between both feet and both hand
    PSD_ratio = transpose(PSD_both_feet./PSD_both_hand);

    figure()
    imagesc(frequencies.Frequencies(band),[1:size(psd.psdt,3)], 10*log10(PSD_ratio))
    c = colorbar;
%     c = colorbar('Ticks',[-50,0,50],...
%                  'TickLabels',{'BF psd << BH psd','BF psd = BH psd','BF psd >> BH psd'});
    c.Limits = [-50,50];
    c.Label.String = 'PSD ratio between both feet and both hand events';
    xlabel('frequency')
    ylabel('channel number')
    set(gca, 'XTick', frequencies.Frequencies(band))
    set(gca, 'YTick', [1:size(psd.psdt,3)])
    title([psd_type{psd_num} ' PSD visualization in the ' band_type{band_num} ' for all channels'])

end
