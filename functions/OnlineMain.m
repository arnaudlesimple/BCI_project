
     %% Data extraction:
     
      clear;
 close all;
    %1) For each gdf file, import data and events
    %Dimension of the data: [NSamples x NChannels]
    addpath(genpath('../projects_common_material/eeglab_current'))
    addpath('C:\Users\Mathieu\Documents\EPFL\Master\MA2\BCI\BCI_project\project2-data-example')
    Map=load('../projects_common_material/channel_location_16_10-20_mi.mat')

    addpath(genpath('../projects_common_material/biosig'));
    filename1 = 'project2-data-example/anonymous.20170613.161402.offline.mi.mi_bhbf.gdf';
    filename2 = 'project2-data-example/anonymous.20170613.162331.offline.mi.mi_bhbf.gdf';
    filename3 = 'project2-data-example/anonymous.20170613.162934.offline.mi.mi_bhbf.gdf'; 
    

    
    [rec1, info1] = sload(filename1);
    [rec2, info2] = sload(filename2);
    [rec3, info3] = sload(filename3);

    % put data from 3 files into one matrix
    AllRec = [rec1;rec2;rec3];
    sizeRec = [length(rec1),length(rec2),length(rec3)];
    Reclabel = [ones(sizeRec(1),1); ones(sizeRec(2),1)*2; ones(sizeRec(3),1)*3];

    AllRecInit = [Reclabel, AllRec(:,[1:16])]; %without reference
    % s: rows correspond to the samples
    %    first column is the filename label
    %    columns 2 to 17 -> 16 channels
    %    column 18 -> reference electrode %je l'enlève

    % Update position data
    info2.EVENT.POS = info2.EVENT.POS + length(rec1);
    info3.EVENT.POS = info3.EVENT.POS + length(rec2) + length(rec1);
    

%%
    load('Classifier\Final classifier.mat')
    load('Classifier\Final classifier Param.mat') %Alpha | Lim
    load('Classifier\Final classifier Features.mat')
    
    
    EEG=AllRec(1:512,1:16);



    test=Online(EEG,classifier,Features_index);
