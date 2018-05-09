function [] = Get_Separability_Plot(Epoch_both_hands, Epoch_both_feet)
% among each frequency and among each channel:
% find the separability distance between hands and feet event

%initialization of the separability matrix (x= freq, y=channel and values
%in boxes refers to separability coeff) [#frequencies x #channels]

%2 boucles for : qui plot psd de both_hands et both_feets over time pour 1 channel et 1 freq
    %i : par channel
    %j : par frequence
    %extraction info i, j from Epoch_both_feet: [time_window x psd]
    %extraction info i, j from Epoch_both_hands [time_window x psd]

    % calcul separability entre psd epoch_both_hands et
    % epoch_both_feet PAR channel PAR freq
    % --> utiliser fonction feature_score (elise)
    %separability_matrix(i,j) = features_score(epoch_both_hands, epoch_both_feet)
    
    


end

