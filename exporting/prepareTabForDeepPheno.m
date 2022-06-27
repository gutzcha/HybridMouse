function [...
    AllClustersForFile,...
    AudioFullName_tSNE_ClusterAnalysis,...
    AllClustersForFileSyllables,...
    AllClustersDurationForFile,...
    TriggerTime,...
    LablesMap] = prepareTabForDeepPheno(tab,saveOptions,LablesMap_new)
%Convert table to cell data variables
%{
Inputs:
        tab - an roi Table or audioClip object
        saveOptions - optional argumant when you want to save the
        variables, a structs with the fields:
            * saveFlag (true/false)
            * fileName (string/char)
            * folder (string/char)
            * saveObj
            * saveExcel
        LablesMap_new (optional) - a labels map, if not inputed, then load
        from default
Outputps:
        AllClustersForFile:
           - cell array with onset of syllable fragmets. each syllable is
           split to 6 ms segments.
        AudioFullName_tSNE_ClusterAnalysis:
           - Path to the audio files
        AllClustersForFileSyllables:
          - cell array with onset of sellables, each syllable with it's own
          length
        AllClustersDurationForFile
         - cell array with duration of each syllable in sec.

If there are syllables from multiple sources, this a cells of cells
%}

%Check input
if isa(tab, 'audioClip')
    obj = tab;
    typeList = tab.allROIs;
    TriggerTime = tab.triggerTime;
    tab = tab.roiTable;
elseif isa(tab,'table')
    obj = [];
    typeList = audioClip.roiTypeList;
    TriggerTime = 0;
else
    error('Invalid input')
end

%Subtract trigger time from tabel

if ~exist('saveOptions','var')||isempty(saveOptions)
    saveOptions.saveFlag = false;
end

if saveOptions.saveFlag == true
    if ~isfield(saveOptions,'fileName')||isempty(saveOptions.fileName)
        if ~isempty(obj)
            saveOptions.fileName = [obj.name,'_data'];
        else
            saveOptions.fileName = 'data';
        end
    end
    if  ~isfield(saveOptions,'folder')||isempty(saveOptions.folder)
        saveOptions.folder = pwd;
    end
    
    if  ~isfield(saveOptions,'saveObj')||isempty(saveOptions.saveObj)||isempty(obj)
        saveOptions.saveObj = false;
    end
    if  ~isfield(saveOptions,'saveExcel')||isempty(saveOptions.saveExcel)||isempty(obj)
        saveOptions.saveExcel = false;
    end
    

    
    
end

if  ~isfield(saveOptions,'subtractTriggerTime')||isempty(saveOptions.subtractTriggerTime)
    saveOptions.subtractTriggerTime = true;
end

%Subtract trigger time
if saveOptions.subtractTriggerTime
    tab.TimeStart = tab.TimeStart-obj.triggerTime;
    tab.TimeEnd = tab.TimeEnd-obj.triggerTime;
end


%Get list of unique paths, each path represents a new file.
AudioFullName_tSNE_ClusterAnalysis = unique(tab.SourcePath); %This is ready
numFiles = size(AudioFullName_tSNE_ClusterAnalysis,1);
numTypes = size(typeList,1);

%Initialize other outputs
AllClustersForFile = cell(numFiles,numTypes);
AllClustersForFileSyllables = cell(numFiles,numTypes);
AllClustersDurationForFile = cell(numFiles,numTypes);

%AllClustersForFile
%Split each syllable into 6 ms fragments
segLenOne = 6e-3; %each segment is 6 ms
%We must seperate between the lines by the file name>cluster(lable)

%Crearte containter - directory between labels and index
LablesMap = (typeList.Name)';
%Function to segment syllables syllables
segmentSyllables= @(s,d,l) s:l:(s+d);

%First seperate into different files #TO DO
for ifile = 1:numFiles
    for itype = 1:numTypes
        thisType = typeList{itype,1};
        inds = tab.Label == thisType;
        start = tab{inds,"TimeStart"};
        duration = tab{inds,"Duration"};
        AllClustersForFileSyllables{ifile,itype} = start;
        AllClustersDurationForFile{ifile,itype} = duration;
        
        % segmentSyllables
        segLen = ones(size(start)).*segLenOne;
        f = arrayfun (segmentSyllables,start,duration,segLen,'UniformOutput',false);
        AllClustersForFile{ifile,itype} =  [f{:,:}]';
        
    end
end

% Get dafult lables map
if ~exist('LablesMap_new','var')||isempty(LablesMap_new)
    path_labelsmap = fullfile('exporting','Labels_Map.csv');
    LablesMap_new = readtable(path_labelsmap);
    LablesMap_new = LablesMap_new.LABEL';
end

[AllClustersForFileSyllables] = sort_labels_in_cellarrays(AllClustersForFileSyllables,LablesMap,LablesMap_new);
[AllClustersDurationForFile] = sort_labels_in_cellarrays(AllClustersDurationForFile,LablesMap,LablesMap_new);
[AllClustersForFile,LablesMap] = sort_labels_in_cellarrays(AllClustersForFile,LablesMap,LablesMap_new);


%Save files
if saveOptions.saveFlag
    saveFullPath = fullfile(saveOptions.folder,saveOptions.fileName);
    if ~saveOptions.saveObj
        save(saveFullPath,'AllClustersForFile','AudioFullName_tSNE_ClusterAnalysis','AllClustersForFileSyllables','AllClustersDurationForFile','TriggerTime','LablesMap','-v7.3')
    else
        save(saveFullPath,'AllClustersForFile','AudioFullName_tSNE_ClusterAnalysis','AllClustersForFileSyllables','AllClustersDurationForFile','TriggerTime','LablesMap','obj','-v7.3')
    end
    P
    if saveOptions.saveExcel
        %Option 1
        %       obj.saveExcle
        %Option 2
        writetable(tab,[saveFullPath,'.xlsx'])
    end
end

function [cellarray_new,LablesMap_new] = sort_labels_in_cellarrays(cellarray_in,LablesMap_old,LablesMap_new)
%Reorder cell array according to a new labels map

%find differences between the two sets
[C,~] = setdiff(LablesMap_old,LablesMap_new);
LablesMap_new = [LablesMap_new, C];

total_num_tags = numel(LablesMap_new);
out_size = [size(cellarray_in,1),total_num_tags];

cellarray_new = cell(out_size);
inds = cellfun((@(x) find(strcmp(x,LablesMap_old))),LablesMap_new,'UniformOutput',false);
for i1 = 1:numel(inds)
    ind = inds{i1};
    if ~isempty(ind)    
        cellarray_new(:,i1) = cellarray_in(:,ind);
    end
end

