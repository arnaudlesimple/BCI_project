clc
clear all
close all
% faire p-e une fonction
%% Extract data
data = load('struct_name')
windows = [] ;
frequencies = [];
electrodes = [];

% struct pour chaque electrode
for i = 1:16
 % extract freq et time par electrode   
end

%% Frequencies VS Baseline

% rapport frequency/baseline

%%
% mettre xlabel, ylabel, nom de la scale bar à ERD / ERS
figure()
subplot(4,4,1)
imagesc(elec1.freq_time)
title('Electrode 1')
colorbar

subplot(4,4,2)
imagesc(elec2.freq_time)
title('Electrode 2')

subplot(4,4,3)
imagesc(elec3.freq_time)
title('Electrode 3')% line plot


subplot(4,4,4)
imagesc(elec4.freq_time)
title('Electrode 4')

subplot(4,4,5)
imagesc(elec5.freq_time)
title('Electrode 5')

subplot(4,4,6)
imagesc(elec6.freq_time)
title('Electrode 6')

subplot(4,4,7)
imagesc(elec7.freq_time)
title('Electrode 7')

subplot(4,4,8)
imagesc(elec8.freq_time)
title('Electrode 8')

subplot(4,4,9)
imagesc(elec9.freq_time)
title('Electrode 9')

subplot(4,4,10)
imagesc(elec10.freq_time)
title('Electrode 10')

subplot(4,4,11)
imagesc(elec11.freq_time)
title('Electrode 11')

subplot(4,4,12)
imagesc(elec12.freq_time)
title('Electrode 12')

subplot(4,4,13)
imagesc(elec13.freq_time)
title('Electrode 13')

subplot(4,4,14)
imagesc(elec14.freq_time)
title('Electrode 14')

subplot(4,4,15)
imagesc(elec15.freq_time)
title('Electrode 15')

subplot(4,4,16)
imagesc(elec16.freq_time)
title('Electrode 16')


%% Plotting data


elec1.time







