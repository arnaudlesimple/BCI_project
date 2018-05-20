function [PosteriorProb] = Online(EEG,classifier,features,decision771,alpha) %alpha

% on recoit de eeg 1S de recordement à 512Hz toute les secondes, ensuite on
% doit les filtrer en frequence et spatiallement, extraire PSD et prédire
% output.

tic

        %CAR filtration
        moyenne=mean(EEG,2);
        EEG =[EEG-moyenne];

    % on garde que les elec qui sont dans nos features
    n_electrode =unique(features(:,2));
    EEG=EEG(:,n_electrode);
        
    features_extracted=[];
    
    for i = 1:size(n_electrode,1)
    
      frequ=features(features(:,2)==n_electrode(i));
      %car il faut au moins deux elements..
      frequ=[frequ',frequ'];
      
      %devrait avoir 16 mais slmnt6
      % on doit faire la psdt de tout car c'est ca la window!!
      
      %[ddd, w, t,psdt1] = spectrogram(EEG(:,i), frequ ,512);
      
      pxx = pwelch(EEG(:,i),[],[],frequ,512) ;

      features_extracted=[features_extracted,pxx(1:end/2)]; %window x frequ
      
    end
    
    score=zeros(1,2);
    
    [predicted_label,score] = predict(classifier, features_extracted); %score [771 773]
    
    
       

           decision771=score(1)*alpha+decision771*(1-alpha);
           
           PosteriorProb=[1-decision771,decision771];

%                 if decision771 > limit771
% 
%                     final_classification=771;
%                 else if decision771(l,m) < limit773
%                       final_classification=773;  
%                     end
%                 end
        
    
    
toc

end

