function [details] = getRecodringDetails(file,path,str)
%[details] = getRecodringDetails(file,path,str)
%file - file name with extention
%path - file path without file name
%str - file name without extention

%   Split the name and set the details


str = lower(str);

%Get rat number
ratNum = str(strfind(str,'rat')+3);

%Get day number
dayNum = str(strfind(str,'day')+3);

%Get paradigm
if contains(str,'chamber')
    paradigm = 'Chamber';
elseif contains(str,'free')
    paradigm = 'Free';
else
    paradigm='';
end

%Get recordingMode (main, mini, sniffing)
if contains(str,'main')
    recordingMode = 'main';
elseif contains(str,'mini')
    recordingMode = 'mini';
elseif contains(str,'sniffing')
        recordingMode = 'sniffing';
else
    recordingMode='';
end

%Get stimulus sex (male, female, mirror)
if contains(str,'female')
    stimulusSex = 'female';
elseif contains(str,'mirror')
        stimulusSex = 'mirror'; %Only one case
else
    stimulusSex='male'; %Default
end

%Get stimulus stimulusAge (Adult (default), juvenile or juv)
if contains(str,'Adult')
    stimulusAge = 'Adult';
elseif contains(str,'juvenile')
    stimulusAge = 'juvenile';
elseif contains(str,'juv')
        stimulusAge = 'juvenile';
else
    stimulusAge='Adult'; %default
end

%Get name of analyst (roee , shira , rotem)
if contains(str,'roee')
    nameOfAnalyst = 'roee';
elseif contains(str,'shira')
    nameOfAnalyst = 'shira';
elseif contains(str,'rotem')
        nameOfAnalyst = 'rotem';
else
    nameOfAnalyst=''; %default
end


details.ratNum = str2double(ratNum);
details.paradigm = paradigm;
details.dayNum = str2double(dayNum);
details.fileName = file;
details.filePath = path;
details.recordingMode = recordingMode;
details.stimulusSex = stimulusSex;
details.stimulusAge = stimulusAge;
details.nameOfAnalyst = nameOfAnalyst;

