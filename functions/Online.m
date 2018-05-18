function [final_classification] = Online(EEG,classifier,features) %alpha

% on recoit de eeg 1S de recordement à 512Hz toute les secondes, ensuite on
% doit les filtrer en frequence et spatiallement, extraire PSD et prédire
% output.

tic


    
    % on garde que les elec qui sont dans nos features
    n_electrode =unique(features(:,2));
    EEG=EEG(:,n_electrode);
        
    features_extracted=[];
    for i = size(n_electrode)
    
      frequ=features(features(:,2)==n_electrode(i));
      %devrait avoir 16 mais slmnt6
      
      [ddd, w, t,psdt1] = spectrogram(EEG(:,i),1 * 512, 0.9375 * 512, frequ ,512);

      features_extracted=[features_extracted,transpose(psdt1)]; %window x frequ
    end
    
    [predicted_label,score] = predict(classifier, features_extracted); %score [771 773]
    

    decision771=0;
    
        for m=1:length(score)

           decision771=score(m)*alpha+decision771*(1-alpha);

                if decision771 > limit771

                    final_classification(m)=771;
                else if decision771(l,m) < limit773
                      final_classification(m)=773;  
                    end
                end
        end
    
    
toc

end

