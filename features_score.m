function [ separation_distance ] = features_score(dataset1, dataset2) 

separation_distance = (mean(dataset1)-mean(dataset2))/sqrt(std(dataset1)+std(dataset2));

end

