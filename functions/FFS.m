function [sel,historyCrit,historyIn] = FFS(selected_data,Action,features,alpha,limit771,limit773,Frequencies,Classifier)

% %% Data aquirement
% clear
% close all
% 
% 
% Action = load('SPD/Event Window.mat');
% Action = Action.Event_window;
% Event=load('SPD/WindowLabel.mat');%WindowLabel
% Event=Event.labelAction;
% 
%  %Attention alpha doit etre penser en terme de practicit� car sinon on a pas le temps
% psd_small_laplacian = load('SPD/SPD with SmallLaplacian Spatial filtre.mat');
% psd_large_laplacian = load('SPD/SPD with LargeLaplacian Spatial filtre.mat');
% psd_CAR_filter = load('SPD/SPD with CAR Spatial filtre.mat');
% psd_no_spatial_filter = load('SPD/SPD with NO Spatial filtre');
% 
% selected_data=psd_CAR_filter.psdt;
% 
% window_frequency = 16;
% frequencies = load('SPD/Frequences.mat');
% load('SPD/Frequences.mat');
% 
% mu_band = 3:6;
% beta_band = 7:18;
% mu_beta_band = 3:18;
% all_band=1:23;
% 
% band = {mu_band, beta_band};
% band_selected = all_band;
% 
% Classifier=[1,2,3,4] %["linear","diaglinear","quadratic","diagquadratic"];
% %%
% %discrimancy = GetDiscrimancyMap(selected_data, band_selected, window_frequency, frequencies);
% %%
% features= [24 9; 12 2; 12 3; 12 5; 12 7;12 8;12 11]; %[frequ x channel]
% 
% alpha=0.1;
% limit771=0.7;
% limit773=0.3;
% 
% 
% 
% %%


maxfeat=size(features,1); %choisit 

   list_features=zeros(1,maxfeat);%1 si choisit 0 sinon
   mieux=0;
   error=zeros(1,maxfeat); %rempli d'erreur
   error(1,1)=1;
   feat_ajoute=zeros(1,maxfeat);
   all_list=zeros(maxfeat,maxfeat);
   
   %% on a fait la premi�re partie youhou!!
    for N_chosen = 1:maxfeat
    
        if mieux==0
         N_run=N_chosen
         
%    selected_inside_index=zeros(1,N_inner);
%    selected_inside_error=zeros(1,N_inner);
%    selected_inside_feature=zeros(1,N_inner);
        
    %loop inner --> TRAIN + VALIDATION
   % for inner_fold = 1:N_inner%loop for Ninner fold inner CV (train+validation error)
        
%         %separate train and validation set
%         train_index = find(partition_inner.training(inner_fold)); %all train indexes
%         validation_index = find(partition_inner.test(inner_fold)); %all validation indexes      
%         train_features = features(train_index,:); %train of inner fold
%         train_labels = labels(train_index); %train of inner fold
%         validation_features = features(validation_index,:); %dont touch until validation error: test of inner fold
%         validation_labels = labels(validation_index); %dont touch until validation error: test of inner fold$
        
        %selected_inside=zeros(1,1);
         
        %je test seulement les features  choisit)
        chosen_features=find(list_features);%chosen features sont ceux qui ont pas des z�ros -> ceux qui on djja �t� choisit
        non_chosen_feat=find(list_features==0); %index/num�ro des features non choissit sur total
        
        inside_error=ones(size(non_chosen_feat));%je vais y mettre les erreurs de chaque ajout
            
        for trie=1:size(non_chosen_feat,2)
            feat_to_try=[chosen_features,non_chosen_feat(:,trie)]; %j'y rajoute un qui n'�tait pas choisit
           [mean_window,std_window,inside_error(trie),std_trial]=GetClassifierAccuracy(selected_data,Action,features(feat_to_try,:),alpha,limit771,limit773,Frequencies,Classifier);
        end
    
        %je choisi le meilleur features a ajouter dans le folder
         selected_inside_error=min(inside_error) %erreur min avec ajout
         index=find(inside_error==min(inside_error)) % si plusieurs fois meme valeur minimumm
         selected_inside_index=index %index du feature dans inside_error = index dans non_selected!
         selected_inside_feature=non_chosen_feat(index) %le num�ro du feature choisis pour chaque fold  
         
         error(N_chosen+1)=selected_inside_error %te donne l'erreur minimum entre tout ajout

   trans=selected_inside_feature;%find(selected_inside_error==min(selected_inside_error,[],2)); %te donne le num�ro du folder ou il y a la min error
   best_index=trans(1,1); %si plusieur min equivalent
   best_feat=best_index; %num�ro du features a rajout� dans  selected
  
  mieux= error(N_chosen+1)>=error(N_chosen) %dit si le nouveau set est meilleur
  
  if mieux==0
  
 
  list_features(best_feat)=1;
  all_list(N_chosen,:)=list_features;
  feat_ajoute(1,N_chosen)=best_feat;
  
  end
  
  if mieux ==1
      Erreur_minimum=error(N_chosen)
      Set_features=feat_ajoute(find(feat_ajoute)) %enleve les z�ros mis lors de l'initialisation
  end
      
        end        
   end       
           
historyCrit=error(2:N_run);
sel=list_features;
historyIn=all_list(1:(N_run-1),:);

end



     
 
