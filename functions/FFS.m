function [selected_features,historyCrit,historyIn] = FFS(selected_data,Action,features,alpha,limit771,limit773,Frequencies,Classifier)




   maxfeat=size(features,1); %choisit 

   list_features=zeros(1,maxfeat);%1 si choisit 0 sinon
   mieux=0;
   error=zeros(1,maxfeat); %rempli d'erreur
   error(1,1)=1;
   feat_ajoute=zeros(1,maxfeat);
   all_list=zeros(maxfeat,maxfeat); %dans le temps
   
   %% on a fait la premi�re partie youhou!!
    for N_chosen = 1:maxfeat
    
        if mieux==0
         N_run=N_chosen
          
        %je test seulement les features  choisit)
        chosen_features=find(list_features);%chosen features sont ceux qui ont pas des z�ros -> ceux qui on djja �t� choisit
        non_chosen_feat=find(list_features==0); %index/num�ro des features non choissit sur total
        
        inside_error=ones(size(non_chosen_feat));%je vais y mettre les erreurs de chaque ajout
            
        for trie=1:size(non_chosen_feat,2)
            feat_to_try=[chosen_features,non_chosen_feat(:,trie)]; %j'y rajoute un qui n'�tait pas choisit
           [inside_error(trie),std_window,mean_trial,std_trial,mean_Confusion]=GetClassifierAccuracy(selected_data,Action,features(feat_to_try,:),alpha,limit771,limit773,Frequencies,Classifier);
        end
    
        %je choisi le meilleur features a ajouter dans le folder en
        %regardant SINGLE SAMPLE ACCURACY
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
           
historyCrit=error(1:N_run); %car error 1=1;
selected_features=list_features;
historyIn=all_list(1:N_run,:);

end



     
 
