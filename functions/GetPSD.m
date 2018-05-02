function GetPSD(beta_band, mu_band, CAR, SmallLaplacian, LargeLaplacian, DoLabel)

    %% Parameter

    if CAR+SmallLaplacian+LargeLaplacian >1 | beta_band+mu_band >1
        A='Error'
        return;
    end
    %% Data extraction:
    %1) For each gdf file, import data and events
    %Dimension of the data: [NSamples x NChannels]
    addpath(genpath('../projects_common_material/eeglab_current'))
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

    %% Isolate interesting bandpass


    % isolate mu and beta bandpass
    sample_rate = info1.SampleRate;

    if mu_band ==1
    % mu band = [7-14] Hz
    [b,a] = butter(4,[7,14]/sample_rate/2);
    fvtool(b,a)
    BandData = filter(b,a,AllRec(:,2:end)); %on prend pas colonne avec label
    AllRecBand=[Reclabel,BandData];
    end

    if beta_band == 1
    % beta band = [15-35] Hz
    [b,a] = butter(4,[15,35]/sample_rate/2);
    fvtool(b,a)
    BandData = filter(b,a,AllRec(:,2:end));
    AllRecBand=[Reclabel,BandData];
    end

    if beta_band ==0 && mu_band==0
        AllRecBand=AllRecInit;
    end

    %% Spatial filtering:

    % CAR filtering:
    %chaque filtre prend AllRec contenant [label, data], freq peut deja être
    %filtré

    if CAR==1   
        car_mu_data = [];
        used_channel=AllRecBand(:,2:end);
        moyenne=mean(used_channel,2);

        FilteredData =[Reclabel, used_channel-moyenne];

        filtre='CAR';
    end

    % Small Laplacian filtering
    if SmallLaplacian ==1   
        small_laplacian = load('small_laplacian.mat');
        FilteredData = [Reclabel, AllRecBand(:,2:end) * small_laplacian.lap];
        filtre='SmallLaplacian';

    end

    % Large laplacian filter
    if LargeLaplacian ==1
        large_laplacian = load('large_laplacian.mat');
        FilteredData = [Reclabel , AllRecBand(:,2:end) * large_laplacian.large_laplacian.lap];

        filtre='LargeLaplacian';
    end

    if LargeLaplacian ==0 && SmallLaplacian==0 && CAR==0
        FilteredData=AllRecBand;

        filtre='NO';
    end

    %filteredDate= [label,datafiltré+band]
    %% Power spectral density 

    time_window = 1; % in seconds
    Init_freq = 512;
    window_change=0.0625;
    time_overlap = 1 - window_change; % in seconds


    frequ_psdt=1/window_change; %frequ overlap

    % 62.5 ms = 16 samples
    Frequencies=4:2:48;

    %look for each recording?

    temps=size(FilteredData(:,:),1)/Init_freq;
    sample=temps*frequ_psdt-(frequ_psdt-1); %car pour les dernieres on ne peut plus décaler!!
    psdt1=zeros(length(Frequencies),sample,16);
    psdt=zeros(sample,length(Frequencies),16);

    window_start=zeros(sample,1);
    window_end=zeros(sample,1);

    for i=1:sample

    window_start(i)=1+(i-1)*(window_change * Init_freq);% = depart absolu*#window*deltaSample -> c'est donc le numéro du sample pas le temps
    window_end(i)=window_start(i)+(time_window*Init_freq)-1;
    end

    for n_electrode = 1:16

      [ddd, w, t,psdt1(:,:,n_electrode)] = spectrogram(FilteredData(:,n_electrode+1),time_window * Init_freq, time_overlap * Init_freq, Frequencies ,512);

      psdt(:,:,n_electrode)=transpose(psdt1(:,:,n_electrode));
    end

    % psdt: columns = frequencies (there are 257 frequencies)
    %       rows = time: each column represents 1 s but we have to take into account the 62.5 ms overlap
    %       -> psdt(:,1) corresponds to PSD estimation between 0 s and 1 s
    %       -> psdt(:,2) corresponds to PSD estimation between 0.0625 s and 1.0625 s
    %       -> psdt(:,3) corresponds to PSD estimation between 0.1250 s and 1.1250 s
    % That is why the number of columns is equal to
    % number_of_samples_of_the_whole_signal * frequ_window/frequ_init - (frequ_window-1)

    %% Epoching

     if DoLabel==1


    BoomTargetHit = 897;
    BoomTargetMiss = 898;
    BothHand = 773; % right cue
    BothFeet=771;
    Fixation = 786;
    ContinuousFeedback = 781;


    StartPositionsHand = (info1.EVENT.POS(info1.EVENT.TYP == BothHand));
    EndPositionsHand = StartPositionsHand+(info1.EVENT.DUR(info1.EVENT.TYP == BothHand));
    Window_BothHand=[];

    Hand_Cue_window=zeros(length(StartPositionsHand),2);

    for i=1:length(StartPositionsHand)
    DebutHand(i)=find(window_end > StartPositionsHand(i),1,'first'); %donne index donc deja le numero de la window!!
    FinHand(i)=find(window_start < EndPositionsHand(i),1,'last');
    Window_BothHand=[Window_BothHand DebutHand(i):FinHand(i)];

    Hand_Cue_window(i,:)=[DebutHand(i),FinHand(i)];
    end

    %Feet
    StartPositionsFeet = (info1.EVENT.POS(info1.EVENT.TYP == BothFeet));
    EndPositionsFeet = StartPositionsFeet+(info1.EVENT.DUR(info1.EVENT.TYP == BothFeet));
    Window_BothFeet=[];

    Feet_Cue_window=zeros(length(StartPositionsFeet),2);

    for i=1:length(StartPositionsFeet)
    DebutFeet(i)=find(window_end > StartPositionsFeet(i),1,'first'); %donne index donc deja le numero de la window!!
    FinFeet(i)=find(window_start < EndPositionsFeet(i),1,'last');
    Window_BothFeet=[Window_BothFeet DebutFeet(i):FinFeet(i)];
    Feet_Cue_window(i,:)=[DebutFeet(i),FinFeet(i)];
    end

    %BoomTargetMiss = 898;

    StartPositionsTargetMiss = (info1.EVENT.POS(info1.EVENT.TYP == BoomTargetMiss));
    EndPositionsTargetMiss = StartPositionsTargetMiss+(info1.EVENT.DUR(info1.EVENT.TYP == BoomTargetMiss));
    Window_BoomTargetMiss=[];

    for i=1:length(StartPositionsTargetMiss)
    DebutTargetMiss(i)=find(window_end > StartPositionsTargetMiss(i),1,'first'); %donne index donc deja le numero de la window!!
    FinTargetMiss(i)=find(window_start < EndPositionsTargetMiss(i),1,'last');
    Window_BoomTargetMiss=[Window_BoomTargetMiss DebutTargetMiss(i):FinTargetMiss(i)];
    end


    %BoomTargetHit = 897;

    StartPositionsTargetHit = (info1.EVENT.POS(info1.EVENT.TYP == BoomTargetHit));
    EndPositionsTargetHit = StartPositionsTargetHit+(info1.EVENT.DUR(info1.EVENT.TYP == BoomTargetHit));
    Window_BoomTargetHit=[];

    for i=1:length(StartPositionsTargetHit)
    DebutTargetHit(i)=find(window_end > StartPositionsTargetHit(i),1,'first'); %donne index donc deja le numero de la window!!
    FinTargetHit(i)=find(window_start < EndPositionsTargetHit(i),1,'last');
    Window_BoomTargetHit=[Window_BoomTargetHit DebutTargetHit(i):FinTargetHit(i)];
    end


    %%
    %ContinuousFeedback = 781;
    StartPositionsFeedBack = (info1.EVENT.POS(info1.EVENT.TYP == ContinuousFeedback));
    EndPositionsFeedBack = StartPositionsFeedBack+(info1.EVENT.DUR(info1.EVENT.TYP == ContinuousFeedback));
    Window_ContinuousFeedback=[];
    Start_End_FeedBack=zeros(length(StartPositionsFeedBack),2);

    %Fixation = 786
    StartPositionsFix = (info1.EVENT.POS(info1.EVENT.TYP == Fixation));
    EndPositionsFix = StartPositionsFix+(info1.EVENT.DUR(info1.EVENT.TYP == Fixation));
    Window_Fixation=[];
    Start_End_Fixation=zeros(length(StartPositionsFix),2);

    Event_window=zeros(length(StartPositionsFix),5);

    for i=1:length(StartPositionsFeedBack)

    DebutFeedBack(i)=find(window_end > StartPositionsFeedBack(i),1,'first'); %donne index donc deja le numero de la window!!
    FinFeedBack(i)=find(window_start < EndPositionsFeedBack(i),1,'last');
    Window_ContinuousFeedback=[Window_ContinuousFeedback DebutFeedBack(i):FinFeedBack(i)];

    DebutFixation(i)=find(window_end > StartPositionsFix(i),1,'first'); %donne index donc deja le numero de la window!!
    FinFixation(i)=find(window_start < EndPositionsFix(i),1,'last');
    Window_Fixation=[Window_Fixation DebutFixation(i):FinFixation(i)];

    type=info1.EVENT.TYP(find(info1.EVENT.POS==(StartPositionsFeedBack(i)))-1);

    Event_window(i,:)=[type,DebutFixation(i) FinFixation(i) DebutFeedBack(i),FinFeedBack(i)];

    end


    %% faire label

    labelAction=zeros(size(psdt,1),1);
    labelAction(Window_BothHand)=773;
    labelAction(Window_BothFeet)=771;
    labelAction(Window_Fixation) = 786;
    labelAction(Window_ContinuousFeedback) = 781;
    labelAction(Window_BoomTargetHit) = 897;
    labelAction(Window_BoomTargetMiss) = 898;

    save('..\SPD\WindowLabel.mat','labelAction');
    save('..\SPD\Frequences.mat','Frequencies');

     end;

    %%
    name= ['..\SPD\SPD with ',filtre,' Spatial filtre.mat'];
    save(name,'psdt');

    save('..\SPD\Feet_Cue_window','Feet_Cue_window');
    save('..\SPD\Hand_Cue_window','Hand_Cue_window');
    save('..\SPD\Event Window','Event_window');

end