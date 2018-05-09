
%Function calculates error of predicted labels compared to real labels of
%samples given in supervised learning
%IMPORTANT: 
% labels1 === label DATA (true labels)
% labels2 === label predicted (by LDA, QDA, ...)

function [ class_error, classification_error ] = error_labels( labels1, labels2)
diff = labels1-labels2; % diff=1 if sample calssified as false when it is true (false negative), diff=-1 if sample classified as true when it is false(false positive). diff=0 if correcly classified.

FalseNeg = sum(diff == 1); %we classified as 0 but is 1 in reality (false negative) "should be one"
FalsePos = sum(diff == -1);%we classified as 1 but is 0 in reality (false positive) "should be zero"

Onetot=sum(labels1==1); %total of positive
Zerotot=sum(labels1==0);%total of negative

class_error= 0.5*(FalseNeg/Onetot)+0.5*(FalsePos/Zerotot);

isCorrect = sum(diff==0); %true positive and true negative
classification_accuracy = isCorrect/(Zerotot+Onetot); %(correctly classified samples)/(total samples)
classification_error = 1-classification_accuracy;
end

%Mathieu: j'ai check et ca a l'air bien