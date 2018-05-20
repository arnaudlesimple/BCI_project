files = 'Mathieu'
addpath(genpath('projects_common_material/biosig'));
if strcmp(files, 'Elise')
    path = ['data/' files '/20180312/'];
elseif strcmp(files, 'Mathieu')
    path = ['data/' files '/20180312/'];
elseif strcmp(files, 'Arnaud')
    path = ['data/' files '/20180319/'];
elseif strcmp(files, 'Anonymous')
    path = ['data/' files '/'];
end

% Initialization
filenames = dir([path '*offline*']);
[AllRec, sizeRec, Reclabel] = deal([],[],[]);
[info1, info2, info3] = deal([],[],[]);
info.EVENT.TYP = [];

for file_num = 1:length(filenames)
    eval(['filename' mat2str(file_num) '= [''' path '''' ' ''' filenames(file_num).name '''];']);
    eval(['[rec' mat2str(file_num) ', info' mat2str(file_num) '] = sload(filename' mat2str(file_num) ');']);
    
    % Put data from 3 files into one matrix
    % -> Define AllRec, sizeRec and Reclabel
    eval(['AllRec = [AllRec; rec' mat2str(file_num) '];']);
    eval(['sizeRec = [sizeRec, length(rec' mat2str(file_num) ')];']);
    eval(['Reclabel = [Reclabel; ones(sizeRec(' mat2str(file_num) '),1) *' mat2str(file_num) '];']); 
end 

AllRecInit = [Reclabel, AllRec(:,[1:16])]; %without reference
% s: rows correspond to the samples
%    first column is the filename label
%    columns 2 to 17 -> 16 channels
%    column 18 -> reference electrode %je l'enlève

% Update position data
if ~isempty(info2)
    info2.EVENT.POS = info2.EVENT.POS + length(rec1);
end
if ~isempty(info3)
    info3.EVENT.POS = info3.EVENT.POS + length(rec2) + length(rec1);
end

info.EVENT.TYP = [];
info.EVENT.POS = [];
info.EVENT.DUR = [];
for file_num = 1:length(filenames)
    eval(['info.EVENT.TYP = [info.EVENT.TYP; info' mat2str(file_num) '.EVENT.TYP];'])
    eval(['info.EVENT.POS = [info.EVENT.POS; info' mat2str(file_num) '.EVENT.POS];'])
    eval(['info.EVENT.DUR = [info.EVENT.DUR; info' mat2str(file_num) '.EVENT.DUR];'])
end