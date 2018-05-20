function GetTopoplot(psd, freq, window_frequency, frequencies, Event_window, window_label, all_plots, all_action, Name)

    addpath(genpath('projects_common_material/eeglab_current'));
    Map=load('projects_common_material/channel_location_16_10-20_mi.mat');
    addpath(genpath('projects_common_material/biosig'));

    
    if all_plots
        n_loop = 4;
        Name={'Small Laplacian','Large Laplacian','CAR','No filter'};
    else
        n_loop = 1;
    end
    
    Action = Event_window;
    Event = window_label.labelAction;

    for a = 1:n_loop
        if all_plots
            PSD = psd{a}.psdt;
        else
            PSD = psd;
        end
        PSD_dB=10.*log10(PSD+1);

        if all_action

            BaseLine=squeeze(mean(PSD_dB((Event==786),:,:),1));
            Hand_Energy=squeeze(mean(PSD_dB((Event==773),:,:),1));
            Feet_Energy=squeeze(mean(PSD_dB((Event==771),:,:),1));

            %for mu band(8:14)
            mu_hand_ERD=mean((Hand_Energy(freq,:)-BaseLine(freq,:))./BaseLine(freq,:),1);
            mu_foot_ERD=mean((Feet_Energy(freq,:)-BaseLine(freq,:))./BaseLine(freq,:),1);
            
            figure()
            B=topoplot(mu_hand_ERD,Map.chanlocs16);
            colorbar;
            title('Hand topoplot');

            figure()
            B=topoplot(mu_foot_ERD,Map.chanlocs16);
            colorbar;
            title('Feet topoplot');

        else


            HandMove=Action(Action==773,:);
            %FeetMove=Action(Action==771,:);
            FeetMove = Action(find(Action(:,1)==771),:);


            for i=1:length(HandMove)
                BaseLineHand=squeeze(mean(PSD_dB(HandMove(i,2):HandMove(i,3),:,:),1));
                HandEpoch=squeeze(mean(PSD_dB(HandMove(i,4):HandMove(i,5),:,:),1));

                BandHandEpoch=mean(HandEpoch(freq,:),1);
                BandBaseLineHand=mean(BaseLineHand(freq,:),1);

                HandERD(i,:)=100*(BandHandEpoch-BandBaseLineHand)./BandBaseLineHand;
                
            end
            %%
            for i=1:length(FeetMove)
                BaseLineFeet=squeeze(mean(PSD_dB(FeetMove(i,2):FeetMove(i,3),:,:),1));
                FeetEpoch=squeeze(mean(PSD_dB(FeetMove(i,4):FeetMove(i,5),:,:),1));

                BandFeetEpoch=mean(FeetEpoch(freq,:),1);
                BandBaseLineFeet=mean(BaseLineFeet(freq,:),1);

                FeetERD(i,:)=100*(BandFeetEpoch-BandBaseLineFeet)./BandBaseLineFeet;
            end



            subplot(2,4,a)
            topoplot(mean(FeetERD),Map.chanlocs16);
            title(['Feet topoplot with ' Name(a) 'filter' ]);
            subplot(2,4,4+a)
            topoplot(mean(HandERD),Map.chanlocs16);
            title(['Hand topoplot with ' Name(a) 'filter' ]);

        end

    end
    
end