clear;
close all;

addpath(genpath('projects_common_material/eeglab_current'));
Map=load('projects_common_material/channel_location_16_10-20_mi.mat');
addpath(genpath('projects_common_material/biosig'));

PSD1 = load('SPD/SPD with CAR Spatial filtre.mat');
PSD1=PSD1.psdt;
PSD2 = load('SPD/SPD with SmallLaplacian Spatial filtre.mat');
PSD2=PSD2.psdt;
PSD3 = load('SPD/SPD with LargeLaplacian Spatial filtre.mat');
PSD3=PSD3.psdt;

AllPSD=[PSD1,PSD2,PSD3];
Name={'CAR','Small Laplacian','Large Laplacian'};

Action = load('SPD/Event Window.mat');
Action = Action.Event_window;
Event=load('SPD/WindowLabel.mat');%WindowLabel
Event=Event.labelAction;

band={'Beta'};
Frequencies=[7:15]; %[3:6](8-14Hz) =mu band, beta= [7:15] (16-32 Hz)

AllAction=0;
SeparateAction=1;

for a=1:3

PSD=AllPSD(:,(a-1)*23+1:a*23,:);
PSD_dB=10.*log10(PSD+1); %+1 in order to be positive

if AllAction==1
    
    BaseLine=squeeze(mean(PSD_dB((Event==786),:,:),1));
    Hand_Energy=squeeze(mean(PSD_dB((Event==773),:,:),1));
    Feet_Energy=squeeze(mean(PSD_dB((Event==771),:,:),1));

    %for mu band(8:14)
    mu_hand_ERD=mean((Hand_Energy(Frequencies,:)-BaseLine(Frequencies,:))./BaseLine(Frequencies,:),1);
    mu_foot_ERD=mean((Feet_Energy(Frequencies,:)-BaseLine(Frequencies,:))./BaseLine(Frequencies,:),1);
    
    subplot(3,2,2*a-1)
    topoplot(mu_foot_ERD,Map.chanlocs16);
    title(['Feet AVERAGE topoplot  of ' band ' waves with ' Name(a) 'filter' ]);
    c=colorbar;
    c.Label.String = 'ERD/ERS';
    
    subplot(3,2,2*a)
    topoplot(mu_hand_ERD,Map.chanlocs16);
    title(['Hand AVERAGE topoplot of ' band ' waves with ' Name(a) 'filter' ]);
    c= colorbar;
    c.Label.String = 'ERD/ERS';
    

    B=topoplot(mu_hand_ERD,Map.chanlocs16);
    colorbar;
    title('Hand topoplot');

    figure;
    B=topoplot(mu_foot_ERD,Map.chanlocs16);
    colorbar;
    title('Feet topoplot');

end

%%


if SeparateAction==1
    
    HandMove=Action(Action==773,:);
    FeetMove=Action(Action==771,:);
 
    
    for i=1:length(HandMove)
         BaseLineHand=squeeze(mean(PSD_dB(HandMove(i,2):HandMove(i,3),:,:),1));
         HandEpoch=squeeze(mean(PSD_dB(HandMove(i,4):HandMove(i,5),:,:),1));
         
         BandHandEpoch=mean(HandEpoch(Frequencies,:),1);
         BandBaseLineHand=mean(BaseLineHand(Frequencies,:),1);
                  
         HandERD(i,:)=100*(BandHandEpoch-BandBaseLineHand)./BandBaseLineHand;

          
    end
%%
    for i=1:length(FeetMove)
         BaseLineFeet=squeeze(mean(PSD_dB(FeetMove(i,2):FeetMove(i,3),:,:),1));
         FeetEpoch=squeeze(mean(PSD_dB(FeetMove(i,4):FeetMove(i,5),:,:),1));
         
         BandFeetEpoch=mean(FeetEpoch(Frequencies,:),1);
         BandBaseLineFeet=mean(BaseLineFeet(Frequencies,:),1);
                  
         FeetERD(i,:)=100*(BandFeetEpoch-BandBaseLineFeet)./BandBaseLineFeet;
         
    end
    

    
    subplot(3,3,1+(a-1)*3)
    topoplot(mean(FeetERD),Map.chanlocs16);
    
    title(['Feet topoplot  of ' band ' waves with ' Name(a) 'filter' ],'FontSize',6);
    c=colorbar;
    c.Label.String = 'ERD/ERS';
    
    subplot(3,3,2+(a-1)*3)
    topoplot(mean(HandERD),Map.chanlocs16);
    title(['Hand topoplot of' band ' waves with ' Name(a) 'filter' ],'FontSize',6);
    c= colorbar;
    c.Label.String = 'ERD/ERS';
    
     subplot(3,3,3+(a-1)*3)
    topoplot(mean(HandERD)-mean(FeetERD),Map.chanlocs16);
    
    title(['Hand-feet topoplot of ' band ' waves with ' Name(a) 'filter' ],'FontSize',6);
    c= colorbar;
    c.Label.String = 'ERD/ERS';
    
   
   
    

    title(['Feet topoplot with ' Name(a) 'filter' ]);
    subplot(3,2,2*a)
    topoplot(mean(HandERD),Map.chanlocs16);
    title(['Hand topoplot with ' Name(a) 'filter' ]);

end

end
