
     %% Data extraction:
     
      clear;
 close all;
    %1) For each gdf file, import data and events
    %Dimension of the data: [NSamples x NChannels]
    addpath(genpath('../projects_common_material/eeglab_current'))
    addpath('C:\Users\Mathieu\Documents\EPFL\Master\MA2\BCI\BCI_project\project2-data-example')
    Map=load('../projects_common_material/channel_location_16_10-20_mi.mat')
    addpath('C:\Users\Mathieu\Documents\EPFL\Master\MA2\BCI\BCI_project - Copie\data\Mathieu\20180312');

    addpath(genpath('../projects_common_material/biosig'));
    filename1 = 'ai9.20180312.142820.offline.mi.mi_bhbf.gdf';
  

    
    [AllRec, info1] = sload(filename1);

%%
   
      B= load('Classifier\supportMathieu');


          %%
    EEG=AllRec(1:512,1:16);    
    lim=size(AllRec,1)/512;
    curr_decision=ones(20,2)*0.5;
    curr_rawprob=zeros(20,2);
    
    for i=2:20
        [curr_decision(i,:),curr_rawprob(i,:)]=Online(AllRec(512*(i-1):512*i,1:16),curr_decision(i-1,2),B);
    end
    %%
    plot(curr_decision(:,2))
    hold on;
    plot(curr_rawprob(:,2))

    