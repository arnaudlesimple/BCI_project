function [curr_decision,curr_rawprob] = Online(EEG,prev_decision,support)

% on recoit de eeg 1S de recordement � 512Hz toute les secondes, ensuite on
% doit les filtrer en frequence et spatiallement, extraire PSD et pr�dire
% output.

tic

%         %CAR filtration
%         moyenne=mean(EEG,2);
%         EEG =[EEG-moyenne];

%small laplacian filtration
        %small_laplacian = load('small_laplacian.mat');
        EEG = EEG * support.small_laplacian;

    % on garde que les elec qui sont dans nos features
    n_electrode =unique(support.selected_features(:,2));
    EEG=EEG(:,n_electrode);
        
    features_extracted=[];
    
    for i = 1:size(n_electrode,1)
    
      frequ=support.selected_features(support.selected_features(:,2)==n_electrode(i));
      %car il faut au moins deux elements..
      frequ=[frequ',frequ'];
      
      %devrait avoir 16 mais slmnt6
      % on doit faire la psdt de tout car c'est ca la window!!
      
      %[ddd, w, t,psdt1] = spectrogram(EEG(:,i), frequ ,512);
      
      pxx = pwelch(EEG(:,i),[],[],frequ,512) ;

      features_extracted=[features_extracted,pxx(1:end/2)]; %window x frequ
      
    end
    
    score=zeros(1,2);
    
    [predicted_label,score] = predict(support.FINALClassifier, features_extracted); %score [771 773]
    
           alpha=support.param(1);

           new_decision=score(1)*alpha+prev_decision*(1-alpha);
           
           curr_decision=[1-new_decision,new_decision];
           
           curr_rawprob=[score(2),score(1)];

           %[773 771]
    
    
toc

end

alpha:.0.01-0.07