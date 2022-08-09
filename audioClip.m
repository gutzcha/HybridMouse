classdef audioClip < handle
    %AUDIOCLIP represent a segment of audiodata
    %   Detailed explanation goes here
    
    properties
        %         + info      - a struct containig info on the source with fields:
        %               * ratNum   - rat serial number/id
        %               * paradigm - experimental paradigm
        %               * dayNum   - day number of experiment
        info = struct('fileName','',...
            'filePath','',...
            'ratNum',0,...
            'paradigm','',...
            'dayNum',0,...
            'recordingMode','main',...
            'stimulusSex','',...
            'stimulusAge','',...
            'nameOfAnalyst','');
        
        vec = [];               % the audio vector, if this is empty, no plot is made
        times_ = [];
        fs = 2.5e5;             % sampling rate of the audio vector
        window = 512;           % fft window length in bins
        overlap = 128;          % fft overlap length in bins
        fft = 512;              % number of DFT point also bins
        threshold = -inf;       % threshold - minimal value for spectrogram
        clim = [];              % color map limits, auto if empty
        climMode = 'auto'       % if clim is on auto, clim will be disregarded
        stft = [];              % short-time Fourier transform of the input signal
        cyclicalF = [];         % cyclicalF   - cyclical frequencies
        timeVec = [];           % vector of time instants
        psd = [];               % psd vector
        spectrogram_sample_rate = [];
        max_bins_to_plot = 2.5e7; %Maximum number of bins to plot
        
        
        emptyFlag = true;        % empty flag, true if the object is empty
        roiTable
        triggerTime = 0;
        sourceType = 'wav'      %type of souce vector - wav or mat
        vecSource = '';        %source of the vector, file of workspace
        ylims;
        absTime = 0;            %If the clip was spliced from a bigger clip, this  property saves the time
        linkRoiTimesToFs = true  %If true, changing sampling rate will also change time stamps.
        detectionParameters = struct('th',0.2,'smoothingK',4,'label',"USV-Auto",'fs',2.5e5)
        cleanSyllables = [];
        probability_vector = array2timetable(nan,'RowTimes',duration(seconds(0))',"VariableNames",{'Probability'});
        load_audio_vector_flag = false; %Determine if the audio vector will be loaded automatically when loading the file
        audioLen = nan;
        obj_snr =struct('mean_snr',nan,'local_snr',nan,'filtered_psd',nan)
        
        
        
    end
    properties (GetAccess = 'public', SetAccess = 'public')
        addedTypeList = table('Size',[0,3],'VariableNames',{'Name','Color','Comments'},'VariableTypes',{'string','double','string'});
    end
    properties (Constant)
        roiTypeList = cell2table(...
            {"USV",[1,0,1],"";...
            'Low',[1,1,0],"";...
            "AC",[0,1,1],"";...
            "USV-Auto",[0,1,0],"AutoDetected USV";...
            "Noise",[1,0,0],"Noise segments identified as USV";...
            "FN",[0,0,0],"Missed syllables";...
            "stimulusUSV",[0,0,1],"Syllables emitted by the stimulus rat";...
            "subjectUSV",[0.5,0,1],"Syllables emitted by the subject rat";...
            "possibly subjectUSV",[1,1,1],"AutoDetected USV by DeepSqueak"},...
            'VariableNames',{'Name','Color','Comments'})
        
        %         roiTableTemplate = table('Size',[0,8],'VariableNames',{'Label','TimeStart','TimeEnd','FrLow','FrHigh','Duration',...
        %             'SourcePath','Comments'},'VariableTypes',{'string','double','double','double','double','double','string','string'});
        roiTableTemplate = table('Size',[0,9],'VariableNames',{'Label','TimeStart','TimeEnd','FrLow','FrHigh','Duration','Score'...
            'SourcePath','Comments'},...
            'VariableTypes',{'string','double','double','double','double','double','double','string','string'});
    end
    properties(Dependent)
        summaryTable
        audioLen_
        allROIs
        path
        times
        name
        recordingMode    %This can be either main or mini or sniffing, derived from the file info
        timeVecAll       % vector of time
        number_of_detected
        extended_roi_table
        
        
    end
    properties(SetAccess = 'private')
        lastDetectionParameters = struct('th',0.4,'smoothingK',4,'label',"USV-Auto","fs",2.5e5)
    end
    methods
        
        function this = audioClip(S,dataIn)
            
            this.roiTable = audioClip.roiTableTemplate;
            
            if nargin==0
                %AUDIOCLIP Construct an instance of this class
                %   Detailed explanation goes here
                return
            else
                if isstruct(S)
                    this.fs            = S.fs           ;
                    this.roiTable      = S.roiTable     ;
                    this.times_        = S.times_       ;
                    this.triggerTime   = S.triggerTime    ;
                    this.ylims         = S.ylims          ;
                    this.absTime       = S.absTime        ;
                    this.sourceType    = S.sourceType     ;
                    this.clim          = S.clim           ;
                    this.climMode      = S.climMode       ;
                    this.window        = S.window       ;
                    this.fft           = S.fft          ;
                    this.overlap       = S.overlap      ;
                    this.addedTypeList = S.addedTypeList;
                    this.vec           = S.vec            ;
                    this.info          = S.info         ;
                    this.audioLen      = S.audioLen;
                    this.obj_snr       = S.obj_snr;
                    if isfield(S,'load_audio_vector_flag')
                        this.load_audio_vector_flag = S.load_audio_vector_flag;
                    end
                else
                    if isvector(S)
                        this.vec = S;
                        this.audioLen = this.audioLen_;
                    end
                    
                    if exist('dataIn','var')&&~isempty(dataIn)
                        this.info = dataIn;
                    end
                end
                this.emptyFlag = false;
            end
        end
        
        function obj_copy = copyobj(obj,p)
            obj_copy = audioClip;
            if ~exist('p','var')||isempty(p)
                p = properties(obj);
            end
            if ~iscell(p)
                p = {p};
            end
            num_att = length(p);
            
            for ip = 1:num_att
                this_att = p{ip};
                mp = findprop(audioClip,this_att);
                % Dont copy constants and dependant properties
                if mp.Dependent || mp.Constant
                   continue 
                end 
                obj_copy.(this_att) = obj.(this_att);
            end
            
        end
        function set.climMode(this,val)
            %Set mode of clim
            %Clim must be a char or string array with two possible values:
            %'auto' or 'manual'
            if ischar(val)||isstring(val)
                
                switch val
                    case {'auto','manual'}
                        this.climMode = val;
                    otherwise
                        error('clim mode must be auto or manual')
                end
            else
                error('invalid clim mode')
            end
        end
        function set.audioLen(this,val)
            if isempty(val)||val==0
                val = this.audioLen_;
            end
            this.audioLen = val;
        end
        function len = get.audioLen_(this)
            
            if isempty(this.vec)
                len = nan;
            else
                len = numel(this.vec)/this.fs;
            end
        end
        
        function [obj,sh,psd] = drawSpectrogram(obj,axSpect, axP)
            %Draw a spectrogram and a normlized pds
            %   Draw plots and update plots accordingly
            
            %{
Inputs:
    axSpect - handel for axis on which to plot the spectrogram
    axP     - handel for axis on which to plot the power spectral density
    dataIn  - an instanse of audioClip class
Outputs:
    sh        - handel for surf plot of spectrogram
    psd       - handel for the line plot psd
    dataOut   - an instanse of audioClip class, updated
       
    
            %}
            if nargin==1
                % Create figure
                figure1 = figure;
                
                % Create axes
                axSpect = axes('Parent',figure1,...
                    'Position',[0.13 0.3 0.775 0.7]);
                
                % Set the remaining axes properties
                set(axSpect,'XTick',zeros(1,0));
                % Create axes
                axP = axes('Parent',figure1,...
                    'Position',[0.13 0.07 0.775 0.23]);
            end
            
            if isempty(obj.vec)
                error ('Empty audio vector')
            end
            if obj.times(1)<0||obj.times(2)>obj.audioLen
                error('Time range must be between %f and %f',0,obj.audioLen)
            end
            vectorInds = (1+obj.fs*obj.times(1)):(obj.fs*obj.times(2));
            %Perform fft on audio vector
            if ~isempty(obj.vec)
                [s,f,t,ps] = spectrogram(obj.vec(vectorInds),obj.window,obj.overlap,obj.fft,obj.fs,'yaxis','MinThreshold',obj.threshold);
            else
                sh = [];
                psd = [];
                
                return
            end
            
            %Create plot with set threshold and clim
            sh = surf(axSpect,log10(abs(s)),'EdgeColor','none');
            view(axSpect,2);
            if ~isempty(obj.clim)
                axSpect.CLim = obj.clim;
            end
            
            %Create psd plot
            psSum = smooth(sum(ps),0.01);
            psSum = psSum./max(psSum);
            psd = plot(axP,t,psSum);
            
            %Make axes tight
            axis(axSpect,'tight')
            axis(axP,'tight')
            
            %Update output values
            obj.stft = s;
            obj.cyclicalF = f;
            obj.timeVec = t;
            obj.psd = ps;
            obj.spectrogram_sample_rate = numel(t)/(t(end)-t(1));
        end
        
        function tab = get.allROIs(obj)
            added = obj.addedTypeList;
            basicTypes = obj.roiTypeList;
            if ~isempty(added)
                tab = [basicTypes;added];
            else
                tab = basicTypes;
            end
            
        end
        function set.fs(obj,newFs)
            %When changing the sampling rate, the time stamps should also be updated to match the new fs
            
            oldFs = obj.fs;
            obj.fs = newFs;
            
            %If the roi table is empty, no need to update anything
            %If there are any syllables in the roi table, all the time
            %stamps should be updated accordingly to the new sampling rate
            
            roiTab = obj.roiTable;
            if isempty(roiTab)||~obj.linkRoiTimesToFs
                return
            end
            
            roiTab.TimeStart = roiTab.TimeStart.*oldFs/newFs;
            roiTab.TimeEnd = roiTab.TimeEnd.*oldFs/newFs;
            obj.roiTable = roiTab;
            obj.refresh_audioLen()
        end
        function set.linkRoiTimesToFs(obj,flag)
            if islogical(flag)&&isscalar(flag)
                obj.linkRoiTimesToFs = flag;
            else
                warning('input to obj.linkRoiTimesToFs must be logical scalar')
            end
        end
        function ret = get.times(obj)
            maxLen = inf;
            if isempty(obj.times_)
                ret(1,1) = 0;
                ret(1,2) =  min(obj.audioLen,1);
            else
                tempRet = obj.times_;
                ret(1) = max(0,tempRet(1));
                ret(2) = min(tempRet(2),obj.audioLen);
                %Make sure the length doen't go over maxLen seconds
                ret(2) = min(ret(2),ret(1)+maxLen);
            end
        end
        function set.times_(obj,t)
            %             maxDisp = 20;
            if isempty(t)
                return
            end
            
            if isnumeric(t)&&(all(size(t) == [1,2]))
                %                 if t(end)<=obj.audioLen&&t(1)>=0&&diff(t)<=maxDisp %#ok
                %                     obj.times_ = t;
                %                 end
                %             else
                %                 error('Invalid time range')
                obj.times_ = t;
            end
        end
        
        function refresh_audioLen(self)
            self.audioLen = self.audioLen_;
        end
        function switch_to_sniffing(self)
            self.fs = 5000;
            self.window = 10;
            self.overlap = 8;
            self.refresh_audioLen()
            
        end
        function set.roiTable(obj,input)
            
            %I have added score variable to the table late, so make sure
            %that older versions get nan values
            if istable(input)
                if ~contains('Score',input.Properties.VariableNames)
                    input.Score = nan(size(input,1),1);
                end
            elseif iscell(input)
                if size(input,2) == 8
                    input = [input,num2cell(nan(size(input,1),2))] ;
                end
            end
            
            try
                [obj.roiTableTemplate;input]; %#ok
            catch ME
                error('Invalid roi table input')
            end
            obj.roiTable = input;
        end
        
        function ret = isempty(obj)
            try
                ret = obj.emptyFlag ;
            catch
                ret = true;
            end
        end
        function set.ylims(obj,vals)
            if all(size(vals)~=[1,2])||~isnumeric(vals)
                %                 warning('ylims must be a 1x2 vector')
                obj.ylims = [0,obj.fs/2]; %#ok
            else
                obj.ylims = vals;
            end
            
        end
        function ret = get.path(obj)
            ret = fullfile(obj.info.filePath,obj.info.fileName);
        end
        function change_folder(obj,new_folder,verbose)
            if ~exist('verbose','var')
               verbose = false; 
            end
            obj.info.filePath = new_folder;
            if verbose
                disp('The path was changed to:')
                disp(obj.path)
            end
        end
        function set.addedTypeList(obj,input)
            
            if isempty(input)
                if ~isempty(obj)
                    %warning('Added roi list was deleted')
                end
                obj.addedTypeList(:,:) = [];
            elseif istable(input)||iscell(input)
                
                
                %If there is an attempt to add a label in the constant
                %label list, dont add it and do not modify
                
                
                constT = obj.roiTypeList;
                
                newNames = input{:,1};
                
                constNames = constT{:,1};
                [~,TF] = setdiff(newNames,constNames);
                input = input(TF,:);
                if isempty(input)
                    return
                end
                newNames = reshape([input{:,1}],1,[]);
                
                %If the new ori types allready exist in the variable,
                %replace them with new comment and new color
                
                t = obj.addedTypeList;
                
                oldNames = [t{:,1}]';
                
                [~,ia,ib] = union(newNames,oldNames);
                t = [t(ib,:);input(ia,:)];
                obj.addedTypeList = t;
            else
                error('Invalid input')
            end
            
        end
        function ret = get.number_of_detected(obj)
            roitab = obj.roiTable;
            ret = size(roitab,1);
        end
        function set.probability_vector(obj,ret)
            if isempty(ret)
                
                return
            end
            try
                prob_vec = ret.Probability;
            catch
                prob_vec = ret.prob_vec;
            end
            if isfield(ret,'detection_times')
                detection_times = ret.detection_times;
            else
                detection_times = [];
            end
            
            if isempty(detection_times)
                obj.probability_vector = [];
            else
                k = linspace(detection_times(1),detection_times(2),(numel(prob_vec)));
                obj.probability_vector = array2timetable(prob_vec','RowTimes',duration(seconds(k))',"VariableNames",{'Probability'});
            end
        end
        
        function recMode = get.recordingMode(obj)
            recMode = obj.info.recordingMode;
        end
        function objName = get.name(obj)
            thisName = obj.info.fileName;
            if isempty(thisName)
                objName = '';
                return
            end
            [~,objName,~] = fileparts(thisName);
        end
        function exlTable = saveExcle(...
                obj,...
                namestring,...
                folderName,...
                subtractTriggerTime,...
                saveFlag,...
                refreshInfoFlag,...
                tableType,...
                fileType,...
                fileName)
            %{
inputs:
    namestring          - a string to add to the saved file, start the name with +
                           or - to add in front or name or before respectivly
    folderName          - save folder name, if empty, save at the pwd
    subtractTriggerTime - if true, subtract the trigger time from all
                          timestamps
    saveFlag            - if true or empty, save the table to disk
    refreshInfoFlag     - if true, reload details from name, not relevant
    tableType           - chooce which table to save, roi for the regular
                          table, extended for extended table and
                          onlyTimeStamps to save only the time stamps
    fileType           -  file type to save table, default : .exls
    fileName           -  file name, if empty, same as the audio file
                %}
                
                if ~exist('saveFlag','var')||isempty(saveFlag)
                    saveFlag = true;
                end
                
                if ~exist('refreshInfoFlag','var')||isempty(refreshInfoFlag)
                    refreshInfoFlag = true;
                end
                
                
                
                if ~exist('subtractTriggerTime','var')||isempty(subtractTriggerTime)
                    subtractTriggerTime = true;
                end
                
                if ~exist('tableType','var')||isempty(tableType)
                    tableType = 'roi';
                end
                if ~exist('fileType','var')||isempty(fileType)
                    fileType = '.xlsx';
                end
                
                if ~exist('fileName','var')||isempty(fileName)
                    fileName = obj.name;
                end
                
                
                if nargin>1
                    if ~isempty(namestring)
                        if namestring(1)=='+'
                            fileName = [fileName,namestring(2:end)];
                        elseif namestring(1)=='-'
                            fileName = namestring(2:end);
                        else
                            fileName = [fileName,namestring];
                        end
                    end
                end
                
                %Add folder name, if empty, will be saved in present folder
                if exist('folderName','var')&&~isempty(folderName)
                    if ~isfolder(folderName)
                        warning('folder %s doesnt exist',folderName)
                    else
                        fileName = fullfile(folderName,fileName);
                    end
                end
                
                %Prevent overwriting
                fileNameCounter = 2;
                fullfilename = [fileName,fileType];
                while isfile(fullfilename)
                    fullfilename = sprintf('%s_%03d%s',fileName,fileNameCounter,fileType);
                    fileNameCounter = fileNameCounter+1;
                end
                
                %In some cases, the delais are missing, regather info to
                %prevent errors
                if refreshInfoFlag
                    obj.info = getRecodringDetails(obj.info.fileName,...
                        obj.info.filePath,obj.info.fileName);
                end
                
                %Add recording data
                tab = obj.roiTable;
                %Subtract trigger time
                if subtractTriggerTime
                    tab.TimeStart = tab.TimeStart-obj.triggerTime;
                    tab.TimeEnd = tab.TimeEnd-obj.triggerTime;
                end
                switch tableType
                    case 'onlyTimeStamps'
                        temp_tab = table(tab.TimeStart,tab.TimeEnd,'VariableNames',{'time_start','time_end'});
                        tab = temp_tab;
                    case 'roi'
                        tab.ratNum = repmat({obj.info.ratNum},[size(tab,1),1]);
                        tab.recordingDay = repmat({obj.info.dayNum},[size(tab,1),1]);
                        tab.paradigme = repmat({obj.info.paradigm},[size(tab,1),1]);
                        tab.stimulusSex = repmat({obj.info.stimulusSex},[size(tab,1),1]);
                        tab.stimulusAge = repmat({obj.info.stimulusAge},[size(tab,1),1]);
                        tab.recordingMode = repmat({obj.recordingMode},[size(tab,1),1]);
                        tab.nameOfAnalyst = repmat({obj.info.nameOfAnalyst},[size(tab,1),1]);
                        %Move file info to front
                        vars = [9:15];
                        tab = movevars(tab,vars,'Before',1);
                    case 'extended'
                        tab = obj.extended_roi_table;
                    otherwise
                        error('The tableType parameter must be either onlyTimeStamps,roi or extended, but it was %s instead',tableType)
                end
                if saveFlag
                    writetable(tab,fullfilename)
                end
                if nargout>0
                    exlTable = tab;
                end
                
        end
        function tabout = filtered_roi_table(obj,labels, times)
            tabout = obj.roiTable;
            if exist('labels', 'var')&&~isempty(labels)
                inds_to_filter = any(tabout.Label== labels,2);
                tabout = tabout(inds_to_filter,:);
            end
            
            if exist('times', 'var')&&~isempty(times)
                inds_to_filter = all([tabout.TimeStart>= times(1),tabout.TimeEnd<= times(2)],2);
                tabout = tabout(inds_to_filter,:);
            end
            


        end
        function tabout = get.extended_roi_table(obj)
            %Reconstcut table to match other programms
            %{
Variables in table:
            
            *. File_Name - this is important when combining several exl
            files
            *. File_Path
            *.Id - serial number of specific usv
            *. Label - the label given to a specific usv
            *. Score - how sure the model is about this usv
            *. Trigger_Time
            *. Begin_Time - start of this usv
            *. End_Time - end of this usv
            *. USV_Length - duration of usv in sec
            *. Mean_Freq
            *. Low_Freq
            *. High_Freq
            *. Range_Freq
            *. Local_SNR
            *. Signal_amp
            *. Signal_PSD
            %}
            
            tabout = cell2table({'','',zeros(0,1),'',zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,125),zeros(0,125),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),...
                zeros(0,1),zeros(0,1),zeros(0,1)},...
                'VariableNames',{'File_Name','File_Path','Id','Label','Score','Local_SNR','Signal_amp','Signal_PSD','Filtered_PSD','Trigger_Time','Begin_Time','End_Time','USV_Length',...
                'Mean_Freq','Low_Freq','High_Freq','Range_Freq'});
            tabout(1,:) = [];
            
            tabin = obj.roiTable;
            
            sz = size(tabin,1);
            if sz==0
                return
            end
            
            File_Name = string(repmat(obj.info.fileName,[sz,1]));
            File_Path = string(repmat(obj.info.filePath,[sz,1]));
            Id = num2cell((1:sz)');
            Label = tabin.Label;
            Score = num2cell(tabin.Score);
            Trigger_Time = num2cell(repmat(obj.triggerTime,[sz,1]));
            Begin_Time = num2cell(tabin.TimeStart);
            End_Time = num2cell(tabin.TimeEnd);
            USV_Length = num2cell(tabin.Duration);
            Mean_Freq = num2cell(mean([tabin.FrHigh,tabin.FrLow],2));
            Low_Freq = num2cell(tabin.FrLow);
            High_Freq = num2cell(tabin.FrHigh);
            Range_Freq = num2cell(abs(tabin.FrHigh-tabin.FrLow));
            Local_SNR = (obj.obj_snr.local_snr);
            Sig_AMP = (obj.obj_snr.sig_amp);
            Sig_psd = (obj.obj_snr.sig_psd);
            Filtered_psd = obj.obj_snr.filtered_psd;
            if size(Local_SNR,1)<sz
                Local_SNR = nan(sz,1);
                Sig_AMP = nan(sz,1);
                Sig_psd = nan(sz,125);
                Filtered_psd = nan(sz,125);
            end
            
            Local_SNR = num2cell(Local_SNR);
            Sig_AMP = num2cell(Sig_AMP);
            Sig_psd = num2cell(Sig_psd);
            tabout = table(File_Name,File_Path,Id,Label,Score,Local_SNR,Sig_AMP,Sig_psd,Filtered_psd,Trigger_Time,Begin_Time,End_Time,USV_Length,Mean_Freq,Low_Freq,High_Freq,Range_Freq);
        end
        function [tab] = roi2table(obj,roilist)
            
            if ~exist('roilist','var')||isempty(roilist)
                tab = obj.roiTableTemplate;
                return
            end
            
            if isa(roilist,'images.roi.Rectangle')
                %                 if isvalid(roilist)
                %
                %                 end
                Label = string(get(roilist,'Label'));
                Score = nan(size(roilist));
                
                %Clear numbers and spaces
                LabelsNoNum = [];
                for ilabel = 1:size(Label)
                    tempLabels = split(Label(ilabel),' ');
                    LabelsNoNum = [LabelsNoNum;tempLabels(1)]; %#ok
                end
                
                Label = LabelsNoNum;
                
                
                
                positions = get(roilist,'Position');
                if iscell(positions)
                    positions = cell2mat(positions);
                end
            elseif isstruct(roilist)&&isfield(roilist,'Positions')
                Label = roilist.Label;
                positions = roilist.Positions;
                if contains('Score',roilist.tab.Properties.VariableNames)
                    Score = roilist.tab.Score;
                else
                    Score = nan(size(roilist.tab,1));
                end
            else
                tab = obj.roiTableTemplate;
                return
            end
            
            
            TimeStart = positions(:,1);
            TimeEnd = TimeStart + positions(:,3);
            FrLow = positions(:,2);
            FrHigh = FrLow + positions(:,4);
            Duration = TimeEnd - TimeStart;
            
            SourcePath = repmat({obj.path},size(TimeStart));
            Comments = repmat("",size(TimeStart));
            
            %             tab = table(Label,TimeStart,TimeEnd,FrLow,FrHigh,Duration,AudioVector,SourcePath,Comments);
            tab = table(Label,TimeStart,TimeEnd,FrLow,FrHigh,Duration,Score,SourcePath,Comments);
            tab = sortrows(tab,{'TimeStart','Label'});
        end
        function obj = replaceRowsInTable(obj,tabIn,deleteAllFlag)
            if ~isa(tabIn,'table')
                tabIn = obj.roi2table(tabIn);
            end
            if ~exist('deleteAllFlag','var')
                deleteAllFlag = false;
                
            end
            if deleteAllFlag
                lables = unique(obj.roiTable.Label);
            else
                lables = unique(tabIn.Label);
            end
            %If there are number, remove them!
            
            
            obj = obj.deleteRowsFromTableByRange(obj.times,lables);
            obj = obj.addRowsToTable(tabIn);
            %             end
        end
        function obj = removeRowsFromTable(obj,lablesToRemove,rangesToRemove)
            %Inputs : lables - labels of rows to remove
            %         ranges - time ranges of rows to remove
            szLabels = size(lablesToRemove,1);
            szRanges = size(rangesToRemove,1);
            if szLabels~=szRanges
                error('Sizes of lables and ranges must be equel')
            end
            if isempty(obj.roiTable)
                return
            end
            
            lablesInObj = obj.roiTable.Label;
            rangesInObj = obj.roiTable{:,{'TimeStart','TimeEnd'}};
            
            indLabels = ismember(lablesInObj, lablesToRemove); %Find inds of tags
            indStarts = ismember(rangesInObj(:,1), rangesToRemove(:,1));
            indEndingss = ismember(rangesInObj(:,2), rangesToRemove(:,2));
            
            %Remove rows where all identifiers are the same
            indToRemove = all([indLabels,indStarts,indEndingss]);
            obj.roiTable(indToRemove,:) = [];
        end
        
        function obj = addRowsToTable(obj,tabIn)
            %Inputs : tabIn - table with data or roi object
            if isempty(tabIn)
                return
            end
            
            
            
            if isa(tabIn,'images.roi.Rectangle')||isa(tabIn,'matlab.graphics.GraphicsPlaceholder')
                tabIn = obj.roi2table(tabIn);
            end
            
            if istable(tabIn)
                if ~contains('Score',tabIn.Properties.VariableNames)
                    tabIn.Score = nan(size(tabIn,1),1);
                end
            elseif iscell(tabIn)
                if size(tabIn,2) == 8
                    tabIn = [tabIn,num2cell(nan(size(tabIn,1),2))] ;
                end
            end
            
            if isempty(obj.roiTable)
                obj.roiTable = tabIn;
            else
                
                obj.roiTable = [obj.roiTable;tabIn];
                obj.roiTable = sortrows(obj.roiTable,{'TimeStart','Label'});
            end
            
        end
        function obj = clearRoiTable(obj,label)
            
            if nargin==1
                obj.roiTable = audioClip.roiTableTemplate;
            else
                t = obj.roiTable;
                t((t.Label == label),:) = [];
                obj.roiTable=t;
            end
        end
        function obj = clearAllAddedROITypes(obj)
            obj.addedTypeList = [];
        end
        function obj = removeDuplicates(obj,also_overlapping,all_overlapping_flag)
            tab = obj.roiTable;
            str = join([tab.Label,string(tab.TimeStart)]);
            [~,ia,~] = unique(str);
            obj.roiTable = audioClip.roiTableTemplate;
            obj.roiTable = tab(ia,:);
            if ~exist('also_overlapping','var')||isempty(also_overlapping)
                also_overlapping = false;
            end
            
            
            if ~also_overlapping
                return
            end
            
            if ~exist('all_overlapping_flag','var')||isempty(all_overlapping_flag)
                all_overlapping_flag = false;
            end
            
            [labels,~,~] = unique(tab.Label);
            if all_overlapping_flag
                temp_tab = remove_duplicates_helper(tab);
            else
                
                temp_tab = audioClip.roiTableTemplate;
                
                for i1 = 1:size(labels,1)
                    this_labels = labels{i1} ;
                    this_tab = tab(tab.Label==this_labels,:);
                    clean_tab = remove_duplicates_helper(this_tab);
                    temp_tab = [temp_tab;clean_tab]; %#ok
                end
            end
            obj.roiTable = temp_tab;
            
            function tab_clean = remove_duplicates_helper(tab,~)
                %                 if ~exist('dominant_tag','var')
                %                    dominant_tag = '';
                %                 end
                
                %From https://uk.mathworks.com/matlabcentral/answers/366626-overlapping-time-intervals
                intervalsIn = [tab.TimeStart,tab.TimeEnd];
                nrow = size(intervalsIn,1);
                % ----- find union of intervals
                intervalsIn = sort(intervalsIn,2);          % make the pairs always increasing pairs
                [intervalsIn, ind] = sort(intervalsIn(:));  % sorts all the input values, and keepts track of how they moved around
                c = [ones(1,nrow) -ones(1,nrow)]';          % this is a matrix that keeps track of when intervals start (1) and stop (-1)
                c = c(ind);                                 % put in order of occurrence
                csc = cumsum(c);                            %sum up starts (1) and stops (-1) , will be =0 at upper end of new interval(s)
                irit = find(csc==0);                        % find index locations of 0 (ends of intervals)
                
                ilef = [1; irit+1];                         % start of intervals index is at the very start (1) and one after all the other ends
                ilef(end) = [];                             % no new interval starting at the very end
                % spansOut matrix is start and end points of the new intervals, y1,y2
                ilef_tab = floor(ilef/2)+1;
                tab_clean = tab(ilef_tab,:);
                tab_clean.TimeStart = intervalsIn(ilef);
                tab_clean.TimeEnd = intervalsIn(irit);
            end
        end
        function obj = deleteRowsFromTableByRange(obj,rangeToDelete,lablesToDelete)
            
            %There is nothing to delete
            if isempty(lablesToDelete)
                return
            end
            
            %There is nothing to delete
            if isempty(obj.roiTable)
                return;
            end
            
            %Out of range
            if ~isempty(rangeToDelete) %If the range is empty, delete every thing
                if rangeToDelete(1)<0||rangeToDelete(2)>obj.audioLen
                    return
                end
            else
                rangeToDelete = [0 obj.audioLen] ;
            end
            
            tab = obj.roiTable;
            %To delet a row is has to have the correct lable and be in the time range
            indsToDelete =  find(any(tab.Label == lablesToDelete',2)&(tab.TimeStart>=rangeToDelete(1) & tab.TimeEnd<=rangeToDelete(2)));
            
            
            if ~isempty(indsToDelete)
                tab(indsToDelete,:)=[];
            end
            obj.roiTable = tab;
            
        end
        
        function ret = autoDetect(obj,data,net,param)
            
            
            if isfield(param,'startTime')
                startTime = param.startTime;
            else
                startTime = obj.times(1);
            end
            
            if isfield(param,'definedFilter')
                definedFilter = param.definedFilter;
            else
                definedFilter = define_a_filter(11,8);
            end
            
            
            if isfield(param,'displayProgressFlag')
                displayProgressFlag = param.displayProgressFlag;
            else
                displayProgressFlag = true;
            end
            
            
            if isfield(param,'fs')
                fsRate = param.fs;
            else
                fsRate = obj.fs;
            end
            
            if isfield(param,'th')
                th = param.th;
                obj.lastDetectionParameters.th = th;
            else
                th = obj.detectionParameters.th;
            end
            
            if isfield(param,'smoothingK')
                smoothingK = param.smoothingK;
                obj.lastDetectionParameters.smoothingK = smoothingK;
            else
                smoothingK = obj.lastDetectionParameters.smoothingK;
            end
            
            if isfield(param,'label')
                label = param.label;
                %Add to added roi list
                obj.addedTypeList = {label,[0 0 0],''};
            else
                label = 'USV-Auto';
            end
            
            ret = extract_features_from_large_filesWraper(...
                data,startTime,definedFilter,net,displayProgressFlag,fsRate,th,smoothingK,label) ;
            
            roiTableTemp = obj.roi2table(ret);
            obj.replaceRowsInTable(roiTableTemp);
            obj.probability_vector = ret;
            
        end
        function explore(obj)
            %If the is no source data, remmber that this was loaded from
            %workspace
            if isempty(obj.path)
                obj.vecSource = 'workspace';
                obj.sourceType = 'mat';
            end
            HybridMouse(obj,audioClip);
        end
        function copyRoi(obj1,obj2,lablesToCopy,rangeToCopy)
            
            if ~exist('lablesToCopy','var')
                lablesToCopy = [];
            end
            
            if ~exist('rangeToCopy','var')||isempty(rangeToCopy)
                rangeToCopy = [0 obj1.audioLen] ;
            end
            
            trigger1 = obj1.triggerTime;
            trigger2 = obj2.triggerTime;
            timeShift = trigger1-trigger2;
            
            T2 = obj2.roiTable;
            T2.TimeStart = T2.TimeStart+timeShift;
            T2.TimeEnd = T2.TimeEnd+timeShift;
            
            obj1.addedTypeList = obj2.addedTypeList;
            if ~isempty(lablesToCopy) %If the specific labels to by copied were specified
                indsToCopy =  contains(T2.Label,lablesToCopy);
                if sum(indsToCopy)>0
                    T2(~indsToCopy,:) = [];
                else
                    warning('Non of the specified labels were found')
                    return
                end
            end
            
            %Remove any roi that exceedes the length of the audio vector
            minStartTime = rangeToCopy(1);
            maxEndTime = rangeToCopy(2);
            
            T2EndTimes = T2.TimeEnd;
            T2StartTimes = T2.TimeStart;
            indsToKeep = (T2StartTimes>=minStartTime)&(T2EndTimes<=maxEndTime);
            if sum(indsToKeep)>0
                T2(~indsToKeep,:)  = [];
            else
                warning ('No syllables within range')
                return
            end
            
            
            obj1.addRowsToTable(T2);
            
        end
        
        function estimatedTriggerTime = findEstimatedTriggerTime(obj,timeRange)
            %estimatedTriggerTime = findEstimatedTriggerTime(obj)
            %Find the estimated trigger time by finding the trigger que -
            %high volt injection in the first 1% of the audio vector
            if ~exist('timeRange','var')||isempty(timeRange)
                len = obj.audioLen;
                timeRange = [0,len*0.01];
                
            end
            
            [~,ind] = max(obj.vecSegment(obj,timeRange)) ;
            estimatedTriggerTime = (ind/obj.fs)+timeRange(1);
        end
        function [currentVec, time_stamps] = getCurrentVec(obj,timeRange,fsIn)
            if ~exist('timeRange','var')||isempty(timeRange)
                timeRange = obj.times;
            end
            if ~exist('fsIn','var')
                fsIn = obj.fs;
            end
            
            currentVec = audioClip.vecSegment(obj,timeRange,fsIn);
            time_stamps = audioClip.get_time_vec(timeRange,fsIn);
        end
        function [vecSegmentOut, newFs, vecSegment] = playSegment(obj,fromFs,toFs,volume,speed,timeRange,vecIn,fsIn)
            % [newVecOut, newFsOut, vecSegment] = playSegment(obj,fromFs,toFs,volume,speed,timeRange)
            
            if ~exist('fromFs','var')||isempty(fromFs)
                fromFs = 60000;
            end
            
            if ~exist('toFs','var')||isempty(toFs)
                toFs = 5000;
            end
            
            
            
            if ~exist('volume','var')||isempty(volume)
                volume = 1;
            end
            
            if ~exist('speed','var')||isempty(speed)
                speed = 0.5;
            end
            
            if ~exist('timeRange','var')||isempty(timeRange)
                timeRange = obj.times;
            end
            
            if ~exist('vecIn','var')||isempty(vecIn)
                vecSegment = obj.vec(round((1+timeRange(1)*obj.fs):(timeRange(2)*obj.fs)));
            else
                vecSegment = vecIn;
            end
            
            if ~exist('fsIn','var')||isempty(fsIn)
                fsIn = obj.fs;
            end
            
            vecSegment = rescale(vecSegment,-1,1);
            
            n = 1;
            outputFs = 44000;
            fsRatio = fromFs/toFs;
            
            %Reduce pitch while maintaning speed
            newFs = fsIn/fsRatio;
            n = n*fsRatio*speed;
            
            %Reduce sampling rate while maintaning pitch and speed
            if newFs>outputFs
                limitFsRatio = newFs/outputFs;
                newFs = outputFs;
                n = n*limitFsRatio;
            end
            
            vecSegmentOutTemp = stretchAudio(vecSegment,n,'Method','wsola');
            if nargout==0
                soundview(vecSegmentOutTemp.*volume,newFs);
            else
                vecSegmentOut = vecSegmentOutTemp;
            end
        end
        
%         function cellArray = getAllSegments(obj,margins, inds, labelsToKeep,randomizeFlag)
%             cellArray =  getAllSegmentsGeneral(obj.vec,obj.fs,obj.roiTable ,margins, inds, labelsToKeep,randomizeFlag);
%         end

        function cellArray = getAllSegments(obj,labelsToKeep,randomizeFlag)
            cellArray = [];
            if (obj.emptyFlag)||isempty(obj.vec) 
                return
            end
          
          
            
            if ~exist('labelsToKeep','var')||isempty(labelsToKeep)
                labelsToKeep = '';
            end
            
            if ~exist('randomizeFlag','var')||isempty(randomizeFlag)
                randomizeFlag = false;
            end
            
            t = obj.roiTable;
            %Filter labels
            if ~isempty(labelsToKeep)
                t(t.Label~=labelsToKeep,:) = [];
            end
                           
            if isempty(t)
                return
            end
            timeRange = [t.TimeStart,t.TimeEnd];
            max_len = obj.audioLen_;
            bad_times = any(timeRange>max_len,2);
            
            if sum(bad_times)>0
                error('Timestamp is out of range at inds %d must not exceed %f',find(bad_times),max_len)
            end
            
            [cellArray, ~] = obj.vecSegment(obj,timeRange,obj.fs);
            
        end
        
        function [axout,out_mat] = plot_psd_surface(obj,normalize_method)
            if ~exist('normalize_method','var')||isempty(normalize_method)
               normalize_method = 1; 
            end
            
            allseg = get_psd(obj,[],normalize_method);
              f = figure;
            ax = axes(f);
        end
        function [axout,out_mat] = plot_psd_mean(obj,normalize_method)
             if ~exist('normalize_method','var')||isempty(normalize_method)
               normalize_method = 1; 
            end
            
            allseg = get_psd(obj,[],normalize_method);
                  f = figure;
            ax = axes(f);
            mean_vec = mean(allseg,1);
            std_vec = std(allseg,1);
            x = 1:numel(mean_vec);
            
            
            x2 = [x, fliplr(x)];
            inBetween = [mean_vec-std_vec, fliplr(mean_vec+std_vec)];
            fill(ax,x2, inBetween, 'g','FaceAlpha',0.1);
            hold on
            plot(ax,x,mean_vec) 
            
%             ax.YLim = [0 125];
            ax.XLim = [0 125];
            ax.XLabel.String = 'Frequency (kHz)';
            ax.YLabel.String = 'Relative Amplitue (A.U.)';
            ax.FontName = 'Times New Roman';
            ax.FontSize = 12;
            
            
            if nargout>0
               axout = ax; 
            end
            if nargout>1
                out_mat = allseg;
            end
        end
        function allseg = get_psd(obj,allseg,normalize_method,filter_flag)
            %normalize_method : 1- normlize per segment, 2 - normalize all
            %0-dont normlize
            %
            
             if ~exist('allseg','var')||isempty(allseg)
                 %Get all audio segment
                allseg = obj.getAllSegments;
             end
             if isempty(allseg)
                 allseg={};
             return
             end
            
            if ~exist('normalize_method','var')||isempty(normalize_method)
                normalize_method = 1;
            end
            
            if ~exist('filter_flag','var')||isempty(filter_flag)
                filter_flag = false;
            end
            
           
            
            %Extract psd from all segment
            allseg = cellfun(@obj.get_power_vector_wrapper,allseg,'UniformOutput',false);
            
            %Reshape to matrix
            allseg = cell2mat(allseg);
            
            %Normlize
            if normalize_method==1
                allseg = normlize_per_file(allseg);
            elseif normalize_method==2
                allseg = normlize_all(allseg);
            elseif normalize_method==0
                return
            else
                error('METHOD must be 1 or 2 but it was %d instead',normalize_method)
            end
          
            
            function mat = normlize_per_file(mat)
                % freq_range = [20:100];
                freq_range = [1:125];
                max_mat = max(mat(:,freq_range),[],2);
                min_mat = min(mat(:,freq_range),[],2);
                mat = (mat-min_mat)./(max_mat - min_mat);
            end
            
            function mat = normlize_all(mat)
                mat = rescale(mat,0,1);
            end
        end
        
        function edit(app)
            splitterEditorApp(app) ;
        end
        function h = editROIs(obj)
            h = editRoisApp(obj);
        end
        function obj1 = plus(obj1,obj2)
            %The must have the same path!
            if ~strcmpi(obj1.path,obj2.path)
                warning ("The objects you are attempting to merg do not have the same file path")
            end
            
            T2 = obj2.roiTable;
            obj1.addedTypeList = obj2.addedTypeList;
            obj1.addRowsToTable(T2);
            
        end
        function [valid_load, errIndo] = setAudioFile(obj,reload_flag,folder_name)
            valid_load = false;
            errIndo = '';
            if ~exist('reload_flag','var')
                reload_flag = false;
            end
            
            if ~exist('folder_name','var')
                folder_name = [];
            end
            
            if reload_flag
                filename = obj.info.fileName;
                
                if ~isempty(folder_name)
                    pathname = folder_name;
                else
                    pathname = obj.info.filePath;
                end
                
                [~,newvec,newinfo,errIndo] = obj.loadAudioFile(filename,pathname);
            else
                [~,newvec,newinfo,errIndo] = obj.loadAudioFile();
                
            end
            
            if ~isempty(newvec)
                obj.vec= newvec;
                obj.info = newinfo;
                obj.audioLen = obj.audioLen_;
                obj.emptyFlag = false;
                valid_load = true;
            else
                warning('The audio vector was not changed!')
            end
        end
        function clean_signal=filterAudio(obj,param)
            if ~exist('param','var')
                param = [];
            end
            
            
            if isfield(param,'noiseTimes')
                noiseTimes = param.noiseTimes;
            else
                noiseTimes = [0 obj.audioLen*0.01]; %if the noise segment was not inputed, assume that first 1% of audio is noise
            end
            
            if isfield(param,'audioTimes')
                audioTimes = param.audioTimes;
            else
                audioTimes = obj.times; %if the noise segment was not inputed, assume that first 10% of audio is noise
            end
            
            
            if isfield(param,'gamma')
                gamma = param.gamma;
            else
                gamma = [] ;
            end
            if isfield(param,'fs')
                fsIn = param.fs;
            else
                fsIn = obj.fs;
            end
            
            
            [noise,noise_bins] = obj.vecSegment(obj,noiseTimes);
            
            if iscell(noise)
                noise = vertcat(noise{:});
                noise_bins = [(noise_bins{:})];
            end
            
            
            num_noise_bins = numel(noise_bins);
            noise_len_sec = num_noise_bins/fsIn;
            
            [vecSegment,~] = obj.vecSegment(obj,audioTimes);
            %             clean_signal=filter_noise(vecSegment, noise,fsIn, gamma) ;
            [clean_signal,~]=WienerNoiseReduction([noise;vecSegment],fsIn,noise_len_sec);
            clean_signal = clean_signal((num_noise_bins+1):end);
            
        end
        function vec = getProbabilityVector(obj,timeRange)
            % vec = getProbabilityVector(obj,timeRange)
            %Create a probability or logic vector
            if ~exist('timeRange','var')||isempty(timeRange)
                timeRange = [0, obj.audioLen]; %Set time range
            end
            
            tab = obj.roiTable;
            
            fs_ = obj.fs;
            abs_onsets = obj.absTime;  %Onset of call (s)
            syl_onsets = obj.roiTable.TimeStart;    %Onset of syllable within a call (ms)
            syl_offsets = obj.roiTable.TimeEnd;  %Onffset of syllable within a call (ms)
            
            
            
            temp_mat_bins = floor([abs_onsets+syl_onsets/1000,abs_onsets+syl_offsets/1000]*fs_); %Abs syl onset and offset in bins
            
            total_len_bin = ceil(total_len*fs_);
            num_of_syllables = size(temp_mat_bins,1);
            
            truth_vec = false(1,total_len_bin); %Total number of bins is = lenght of output
            
            for i1 =1:num_of_syllables
                inds = 1+temp_mat_bins(i1,1):temp_mat_bins(i1,2);
                if numel(inds)==0
                    inds = temp_mat_bins(i1,1)+1;
                end
                truth_vec(inds)=true  ;
            end
            
        end
        function set.obj_snr(obj,val)
            ret = struct('mean_snr',nan,'local_snr',nan,'sig_amp',nan,'sig_psd',nan(1,125),'filtered_psd',nan(1,125));
            if isempty(val)
                obj.obj_snr = ret;
                return
            else
                if isfield(val,'mean_snr')
                    ret.mean_snr = val.mean_snr;
                end
                if isfield(val,'local_snr')
                    ret.local_snr = val.local_snr;
                end
                  if isfield(val,'filtered_psd')
                    ret.filtered_psd = val.filtered_psd;
                  end
                 if isfield(val,'sig_amp')
                    ret.sig_amp = val.sig_amp;
                 end
                  if isfield(val,'sig_psd')
                    ret.sig_psd = val.sig_psd;
                end
                obj.obj_snr = ret;
                
            end
        end

        function ret = prepareTabForDeepPhenoObj(obj,saveOptions)
            if ~exist('saveOptions','var')
                saveOptions = [];
            end
            [AllClustersForFile,...
                AudioFullName_tSNE_ClusterAnalysis,...
                AllClustersForFileSyllables,...
                AllClustersDurationForFile,TriggerTime,LablesMap] = prepareTabForDeepPheno(obj,saveOptions);
            ret.AllClustersForFile = AllClustersForFile;
            ret.AudioFullName_tSNE_ClusterAnalysis = AudioFullName_tSNE_ClusterAnalysis;
            ret.AllClustersForFileSyllables = AllClustersForFileSyllables;
            ret.AllClustersDurationForFile = AllClustersDurationForFile;
            ret.TriggerTime = TriggerTime;
            ret.LablesMap = LablesMap;
            
        end
        function rename_label(this,old_label,new_label)
            arguments
                this
                old_label (1,:) {mustBeText}
                new_label (1,:)
            end
            t = this.roiTable;
            t = this.roiTable;
            
            inds = t.Label == old_label;
            num_to_replace = sum(inds);
            
            if isempty(new_label)
                t(inds,:) = [];
                fprintf('%d row with %s labels were removed\n',num_to_replace,old_label)
            else
                
                if ischar(old_label)
                    old_label = string(old_label);
                end
                if ischar(new_label)
                    new_label = string(new_label);
                end
                
                
                t.Label(inds) = new_label;
                if num_to_replace
                    fprintf('%d rows with %s labels were replaced by %s labels\n',num_to_replace,old_label,new_label)
                else
                    fprintf('%s label was not found, no rows were replaced\n',old_label)
                end
            end
            this.roiTable = t;
        end
        
        function [ret] = refresh_obj_snr(obj,calculate_sepratly_flag, filter_flag)
            % Calculate local SNR for each label seperatly
            
            % info_table: a table specifing the important labels and the
            % frequency ranges
            ret.mean_snr = nan;
            ret.local_snr = nan;
            ret.filtered_psd = nan;
            ret.sig_amp = nan;
            ret.sig_psd = nan;
            
            if isempty(obj)
                warning('The object is empty')
                return
            end
            if obj.emptyFlag
                warning('The object is empty')
                return
            end
            if isempty(obj.vec)
                warning('No audio vector loaded')
                return
            end
            
            if ~exist('calculate_sepratly_flag','var')||isempty(calculate_sepratly_flag)
                calculate_sepratly_flag = true;
            end
            
            if ~exist('filter_flag','var')||isempty(filter_flag)
                filter_flag = true;
            end
            
          
            
            
            audio_length = obj.audioLen_;
            if isnan(audio_length)
                warning('No audio vector loaded')
                return
            end
            
            tab = obj.roiTable;
            all_labels = unique(tab.Label);
            num_all_labels = length(all_labels);
%             labels_to_omit = ["Noise"];
            labels_to_process = ["Low","subjectUSV","stimulusUSV"];
            ommited_labels = all_labels(~contains(all_labels,labels_to_process));
            all_labels= all_labels(contains(all_labels,labels_to_process));
           
            num_labels = length(all_labels);
            
            variable_names = {'Label','TimeStart','TimeEnd','mean_snr','local_snr','filtered_psd','sig_amp','sig_psd'};
            summary_tab = cell2table(cell(num_all_labels,numel(variable_names)),"VariableNames",variable_names);
            if calculate_sepratly_flag %TODO make it possible to extract the parameters on all labels together
                
            else
                
            end
            
            for ilabel = 1:num_labels
                label = all_labels{ilabel};
                
                %Adjust "low" filter frequency
                tab_temp = tab(tab.Label==label, :);
                if label=="Low"
                    frlow = tab_temp.FrLow;
                    frlow(frlow<2000) = 2000;
                    
                    
                    frhigh = tab_temp.FrHigh;
                    frhigh(frhigh<4000) = 4000;
                    
                    fix_inds = frhigh<=frlow;
                    frhigh(fix_inds) = frlow(fix_inds)+3000;
                    
                    tab_temp.FrLow = frlow;
                    tab_temp.FrHigh = frhigh;
                    
                end
                
                
                ret_temp = get_spectral_data(obj, tab_temp, audio_length, filter_flag);
                summary_tab{ilabel,'Label'} = {label};
                summary_tab{ilabel,'TimeStart'} = {tab_temp.TimeStart};
                summary_tab{ilabel,'TimeEnd'} = {tab_temp.TimeEnd};
                summary_tab{ilabel,'mean_snr'} = {ret_temp.mean_snr};
                summary_tab{ilabel,'local_snr'} = {ret_temp.local_snr};
                summary_tab{ilabel,'filtered_psd'} = {ret_temp.filtered_psd};
                sig_amp_temp = ret_temp.sig_amp;
                if any(isnan(sig_amp_temp))
                    aa=  1
                end
                summary_tab{ilabel,'sig_amp'} = {ret_temp.sig_amp};
                summary_tab{ilabel,'sig_psd'} = {ret_temp.sig_psd};
                
            end
            
            num_labels_to_omit = length(ommited_labels);
            for ilabel = (1+num_labels):(num_labels_to_omit+num_labels)
                label = ommited_labels(ilabel-num_labels);
                tab_temp = tab(tab.Label==label, :);
                summary_tab{ilabel,'Label'} = {label};
                summary_tab{ilabel,'TimeStart'} = {tab_temp.TimeStart};
                summary_tab{ilabel,'TimeEnd'} = {tab_temp.TimeEnd};
                summary_tab{ilabel,'mean_snr'} = {nan};
                summary_tab{ilabel,'local_snr'} = {nan(numel(tab_temp.TimeStart),1)};
                summary_tab{ilabel,'filtered_psd'} = {nan(numel(tab_temp.TimeStart),125)};
                summary_tab{ilabel,'sig_amp'} = {nan(numel(tab_temp.TimeStart),1)};
                summary_tab{ilabel,'sig_psd'} = {nan(numel(tab_temp.TimeStart),125)};
            end
            
            % Reorgnize
            len_table = height(summary_tab);
            labels = [];
            time_start = [];
            time_end = [];
            mean_snr = cell(0,2);
            local_snr = [];
            filtered_psd = double.empty(0,125);
            sig_amp = [];
            sig_psd = double.empty(0,125);
            for itab = 1:len_table
                temp_times_start = summary_tab{itab,'TimeStart'}{:};
                num_syl = length(temp_times_start);
                temp_labels = repmat(summary_tab{itab,'Label'}{:},num_syl,1);
                labels = [labels;string(temp_labels)];
                time_start = [time_start;temp_times_start];
                time_end = [time_end;summary_tab{itab,'TimeEnd'}{:}];
                mean_snr = [mean_snr ; {summary_tab{itab,'Label'}{:},summary_tab{itab,'mean_snr'}{:}}];
                local_snr = [local_snr;summary_tab{itab,'local_snr'}{:}];
                filtered_psd = [filtered_psd;summary_tab{itab,'filtered_psd'}{:}];
                sig_amp = [sig_amp;summary_tab{itab,'sig_amp'}{:}];
                sig_psd = [sig_psd;summary_tab{itab,'sig_psd'}{:}];
                    
            end
            all_data = table(labels,time_start,time_end,local_snr,filtered_psd,sig_amp,sig_psd,'VariableNames',{'Label','TimeStart','TimeEnd','local_snr','filtered_psd','sig_amp','sig_psd'});
            all_data = sortrows(all_data,'TimeStart');
            %Record this
            ret.mean_snr = mean_snr;
            ret.local_snr = all_data.local_snr;
            ret.filtered_psd = all_data.filtered_psd;
            ret.sig_amp = all_data.sig_amp;
            ret.sig_psd = all_data.sig_psd;
            
            
            obj.obj_snr = ret;
            
            function [ret, isi_times] = get_spectral_data(obj,tab,audio_length, filter_flag)
                   %             tab(tab.Label='Noise',:) = [];
            
            if ~exist('filter_flag','var')||isempty(filter_flag)
                filter_flag = false;
            end
            signal_requency_range = [tab.FrLow, tab.FrHigh];
            signal_times = [tab.TimeStart,tab.TimeEnd];
            if isempty(signal_times)
                return
            end
            
            signal_times(signal_times>audio_length) = audio_length;
            
            %validate times
            if any(diff(signal_times,[],2)<=0)
                warning('Something is wrong with the time table')
            end
            
            
            
            isi_times = get_isi_times(signal_times,audio_length);
            try
                audio_vec = audioClip.vecSegment(obj,signal_times);
                % if there is only one syllable, the output is a double array,
                % but this has to be a cell array
                if ~iscell(audio_vec)
                   audio_vec = {audio_vec};
                end
            catch ME
                disp(ME.message)
            end
            noise_vec = audioClip.vecSegment(obj,isi_times);
               
            if filter_flag
               % band filter each signal 
               
               %Remove outliers and determine filter range
               frlow_temp = filloutliers(tab.FrLow','clip','quartiles');
               frhigh_temp = filloutliers(tab.FrHigh','clip','quartiles');
               min_fr = min(frlow_temp);
               max_fr = max(frhigh_temp);
               
               %filter function  
               
               
               filter_func = @(x) bandpass(x,[min_fr,max_fr],obj.fs);
               
               %filter the audio and noise signals before calculating 
               audio_vec = cellfun(filter_func,audio_vec,'UniformOutput',false);
               noise_vec = cellfun(filter_func,noise_vec,'UniformOutput',false);
            end
            %get_amp = @(y) sqrt(mean(y.^2));
            %get_amp = @(y) bandpower(y,2.5e5,[0 2.5e5/2]);
            %get_amp = @(y) rms(y)^2;
            get_amp = @(y) db(periodogram(y,kaiser(length(y),38) ,obj.fs*2-1,obj.fs));
            
            signal_noise_db = (cellfun(get_amp,audio_vec,'UniformOutput',true));
            noise_db =  (cellfun(get_amp,noise_vec,'UniformOutput',true));
            noise_db(isnan(noise_db)) = mean(noise_db,'omitnan'); %Fill missing values
           
            
            ret.sig_amp = signal_noise_db;
            
            ret.local_snr = get_local_snr(signal_noise_db,noise_db);
            
            ret.sig_psd = obj.get_psd(audio_vec,1);
            noise_psd = obj.get_psd(noise_vec,1);
            nan_inds = any(isnan(noise_psd),2);
            if sum(nan_inds)>0
                noise_psd(nan_inds,:) = repmat(mean(noise_psd,1,'omitnan'),[sum(nan_inds),1]); %Fill missing values
            end
            ret.filtered_psd = get_local_psd(ret.sig_psd,noise_psd);
            
            ret.mean_snr =  get_weighted_mean(audio_vec,signal_noise_db,noise_vec,noise_db);

            end
            
            function mean_snr = get_weighted_mean(audio_vec,signal_noise_db,noise_vec,noise_db)
                numels_signal = cellfun(@numel,audio_vec);
                
                % remove files with problems
                empty_ind = numels_signal==0;
                numels_signal(empty_ind,:) = [];
                signal_noise_db(empty_ind) = [];
                
                weight_signal = numels_signal./sum(numels_signal);
                wighted_mean_signal = sum(signal_noise_db.*weight_signal);
                
                numels_noise = cellfun(@numel,noise_vec);
                weight_noise = numels_noise./sum(numels_noise);
                wighted_mean_noise = sum(noise_db.*weight_noise);
                
                
                mean_snr = (wighted_mean_signal-wighted_mean_noise)/wighted_mean_noise;
            end
            function local_snr = get_local_snr(sig,isi)
                sz = size(sig,1);
                local_snr = zeros(sz,1);
                for i1 = 1:sz
                    local_noise = mean(isi(i1:i1+1));
                    local_snr(i1) = (sig(i1) - local_noise)/local_noise;
                end
                
            end
                function local_psd = get_local_psd(sig,isi)
                sz = size(sig,1);
                local_psd = zeros(sz,125);
                for i1 = 1:sz
                    local_noise = mean(isi(i1:i1+1,:));
                    local_psd(i1,:) = (sig(i1,:) - local_noise);
                    if any(isnan(local_psd(i1,:)))
                       a = 1; 
                    end
                end
                
            end
            function isi_times = get_isi_times(signal_times,audio_length)
                sz = size(signal_times,1);
                isi_times = [signal_times(1:sz-1,2),signal_times(2:sz,1)];
                start_seg(2) = max(signal_times(1,1),0);
                start_seg(1) = max(start_seg(2)-20,0);
                if diff(start_seg)<-0
                    start_seg = signal_times(1,:);
                end
                
                end_seg(2) = min(signal_times(end,2)+20,audio_length);
                end_seg(1) = min(signal_times(end,2),audio_length);
                if diff(end_seg)<=0
                    end_seg =signal_times(end,:);
                end
                
                isi_times = [start_seg;isi_times;end_seg];
                
            end
            
     
        end
        function plot_current_vec_with_vocalization(obj, time_range, labels)
            if ~exist('time_range','var')||isempty(time_range)
                time_range = obj.times;
            end
             if ~exist('labels','var')
                labels = [];
             end
             
            [vec_,t] = obj.getCurrentVec(time_range);
            tab = obj.filtered_roi_table(labels,time_range)  ;
            voc_times = [tab.TimeStart,tab.TimeEnd];
            t_voc = nan(size(t));
            for it = 1:size(voc_times,1)
                t_voc((voc_times(it,1)<t) & (t<=voc_times(it,2))) = vec_((voc_times(it,1)<t) & (t<=voc_times(it,2)));
            end
            f = figure;
%             ax = f.CurrentAxes;
            plot( t,vec_)
            hold on
            plot( t,t_voc)
            hold off
            title(obj.name)
            ax = f.CurrentAxes;
            ax.XLabel.String = 'Time (sec)';
%             ax.YLabel.String = ''
        end
        function denoised=get_power_vector_wrapper(obj,signal_data,noise_data,filter_flag)
                
           
                normalization_frequency_range = [40,120];
                if ~exist('filter_flag','var')
                    filter_flag = false;
                end
                
                
                %Get power spectral density of the signal and the noise
                frequency_range = [1,125];
                
                signal_b = get_power_vector(signal_data,frequency_range);
                if filter_flag
                noise_b = get_power_vector(noise_data,frequency_range);
                else
                   noise_b = zeros(size(signal_b)) ;
                end
                
                %Subtract
                denoised = signal_b-noise_b;
                
                
                %Normlize
                min_val = min(denoised(normalization_frequency_range(1):normalization_frequency_range(end)),[],"all");
                max_val = max(denoised(normalization_frequency_range(1):normalization_frequency_range(end)),[],"all");
               
                denoised=(denoised-min_val)/(max_val-min_val);
              
                
            end
              
        function [new_tab,old_tab] = refresh_scores(obj)
            
            t = obj.roiTable;
            old_t = t;
            prob_vec = obj.probability_vector.Probability;
            tab.Onsets = t.TimeStart;
            tab.Offsets = t.TimeEnd;
            fs_new = (numel(prob_vec))/obj.audioLen;
            
            scores = get_scores_from_probabilty_vector(...
                prob_vec,tab,fs_new);
            t.Score = scores;
            obj.roiTable = t;
            disp('Scores were refreshed')
            if nargout>1
                new_tab = t;
            end
            if nargout>2
                old_tab = old_t;
            end
        end
        
        
        function s = saveobj(this)
            s.info          = this.info         ;
            s.vec           = []                ;
            s.times_        = this.times_       ;
            s.fs            = this.fs           ;
            s.window        = this.window       ;
            s.overlap       = this.overlap      ;
            s.fft           = this.fft          ;
            s.threshold     = this.threshold    ;
            s.clim          = this.clim         ;
            s.climMode      = this.climMode     ;
            s.path          = this.path         ;
            s.stft          = []                ;
            s.cyclica       = []                ;
            s.timeVec       = []                ;
            s.psd           = []                ;
            
            s.addedTypeList = this.addedTypeList;
            s.emptyFlag     = this.emptyFlag    ;
            s.roiTable      = this.roiTable     ;
            s.sourceType    = this.sourceType   ;
            s.ylims         = this.ylims        ;
            s.absTime       = this.absTime      ;
            s.triggerTime   = this.triggerTime  ;
            s.load_audio_vector_flag = this.load_audio_vector_flag;
            s.audioLen      = this.audioLen;
            s.obj_snr       = this.obj_snr;
            
        end
        
    end
    methods(Static)
        
        function tab_out = construct_roi_from_timestamps(times, Label)
            tab_out = audioClip.roiTableTemplate;
            if ~exist('times','var')||isempty(times)
               return
            end
            if ~exist('labels','var')||isempty(Label)
               Label =  'USV';
            end
            if istable(times)
               times = table2array(times) ;
            end
            sz = size(times,1);
            sz_labels = size(Label,1);
            if sz_labels>1 && size(label,1)~=sz
               warning ('The number of labels must be 1 or match the height of the times but it was %d instead',sz_labels)
               return
            elseif sz_labels==1
                Label = repmat(string(Label),sz,1);
            end
         
            TimeStart = times(:,1);
            TimeEnd = times(:,2);
            FrLow = nan(sz,1);
            FrHigh = nan(sz,1);
            Duration = TimeEnd-TimeStart;
            Score = nan(sz,1);
            SourcePath = repmat("",sz,1);
            Comments = repmat("",sz,1);
            tab_temp = table(Label,TimeStart,TimeEnd,FrLow,FrHigh,Duration,Score,SourcePath,Comments);
            tab_out = [tab_out;tab_temp];
            
        end
        
        function [vecseg, all_inds] = vecSegment(obj,timeRange,fs)
            vecseg = [];
            all_inds = [];
            if obj.emptyFlag
                
                return
            end
            if isempty(obj.vec)
                return
            end

            if ~exist('timeRange','var')
                timeRange = obj.times;
            end
            
            [numSegments, checkInput] = size(timeRange);
            if checkInput~=2
                error('timeRange must have size of nx2 (has %dx%d)',numSegments, checkInput)
            end
            
            vecseg = cell(numSegments,1);
            all_inds = cell(numSegments,1);
            for iSeg = 1:numSegments
                if isnumeric (obj)
                    if ~exist('fs','var')
                        fs = 2.5e5;
                        warning ('Sampling rate was set to %0.1f',fs);
                    end
                    inds = round((1+fs*timeRange(iSeg,1)):(fs*timeRange(iSeg,2)));
                    
                    vecseg{iSeg,1} = obj(inds);
                    all_inds{iSeg,1} = inds;
                    
                elseif isa(obj,'audioClip')
                    inds = round((1+obj.fs*timeRange(iSeg,1)):(obj.fs*timeRange(iSeg,2)));
                    vecseg{iSeg,1} = obj.vec(inds);
                    all_inds{iSeg,1} = inds;
                else
                    warning ('Invalid input type %s',class(obj))
                end
                
            end
            if numSegments==1
                vecseg = vecseg{:,:};
                all_inds = all_inds{:,:};
            end
        end
        function timevec = get_time_vec(time_range,fs_)
            
            sz = size(time_range,1);
            if sz==1
         
            time_diff = diff(time_range);
            n_points = round(time_diff*fs_)+1;
            timevec = linspace(time_range(1),time_range(2),n_points);
            timevec(end) = [];
            else
                all_vec = cell(sz,1);
                for i1 = 1:sz
                    all_vec{i1} = audioClip.get_time_vec(time_range(i1,:),fs_);

                end
                timevec = all_vec;
            end
        end
        function [file,vec,info,errInfo] = loadAudioFile(file,path)
            vec = [];
            info = [];
            errInfo = '';
            if nargin ==0
                [file,path] = uigetfile(...
                    {'*.mat;*.wav','All Files (*.mat,*.wav)';...
                    '*.wav','Audio File (.*wav)';...
                    '*.mat','MATLAB File (*.mat)'});
            end
            
            
            if file==0
                warning('No file was chocen, returning an empty audioClip file')
                return
            else
                try
                    [~,f,type]=fileparts(file);
                    switch lower(type)
                        case '.wav'
                            vec = audioread(fullfile(path,file)) ;
                        case '.mat'
                            readfile = matfile(fullfile(path,file));
                            classesInVar = whos(readfile);
                            indDouble = find(strcmpi('double',{classesInVar.class}),1,'first');
                            vec  = eval(['readfile.',classesInVar(indDouble).name]);
                    end
                    info = getRecodringDetails(file,path,f);
                catch ME
                    warning('Unable to find the audio file on Matlab path, please call the function again without argumants to reload the file manualy')
                    disp(ME.message)
                    errInfo = ME.message;
                end
            end
        end
        
        function this = loading_options(S)
            answer = questdlg('Unable to find file on path, please choose how to proceed', ...
                'Load audio file', ...
                'Load Manualy','Load without audio vector','Cancel','Cancel');
            
            switch answer
                case 'Load Manualy'
                    [file,vec,info] = audioClip.loadAudioFile();
                    if file==0
                        %No file was selected
                        this = audioClip;
                        return;
                    else
                        %Set new info and vector
                        S.info = info;
                        S.vec = vec;
                        if isstruct(S)
                            this = audioClip(S);
                        elseif isa(S,'audioClip')
                            this = S;
                        end
                        return
                    end
                    
                case 'Load only roiTable'
                    this = S.roiTable;
                case 'Load without audio vector'
                    S.vec = [];
                    this = audioClip(S);
                case 'Cancel'
                    %return empty audioClip
                    this = audioClip;
                    return
            end
        end
        
        
        function this = loadobj(S)
            %Determine if to load audio vector (true) or only roi table (false)
            %             flags = false;
            if ~isfield(S,'load_audio_vector_flag')
                flags = false;
            else
                flags = S.load_audio_vector_flag;
            end
            
            if isstruct(S)
                path = fullfile(S.path);
                if ~isfield(S,'sourceType')
                    S.sourceType = 'wav' ;
                end
                if ~isfield(S,'ylims')
                    S.ylims = [];
                end
                
                if ~isfield(S,'absTime')
                    S.absTime = 0;
                end
                
                if ~isfield(S,'clim')
                    S.clim = [];
                end
                
                if ~isfield(S,'climMode')
                    
                    S.climMode = 'auto';
                end
                if ~isfield(S,'triggerTime')
                    S.triggerTime = 0;
                end
                if ~isfield(S,'audioLen')
                    S.audioLen = [];
                end
                
                %Clear/fix Score from roi table
                roitab = S.roiTable;
                
                if ~contains('Score',roitab.Properties.VariableNames)
                    roitab.Score = nan(size(roitab,1),1);
                else
                    if size(roitab.Score,2) > 1
                        roitab.Score = nan(size(roitab,1),1);
                    end
                end
                
                if ~isfield(S,'obj_snr')
                    S.obj_snr = [];
                end
                
                
                S.roiTable = roitab;
                
                %Make sure file name has right extention
                [~,fileName,ext] = fileparts(path);
                if isempty(ext)
                    ext = S.sourceType ;
                    path = [path,ext];
                else
                    if ~strcmpi(ext(2:end),S.sourceType)
                        warning('File extention (%s) and sourceType (%s) are not the same,\nChaging to extention',ext,S.sourceType)
                        S.sourceType = ext(2:end);
                    end
                end
                
                %Search for the file with full path, if it cant be found,
                %try to look for the file in current folder
                if ~isfile(path)
                    %if the file was not found then search for the file
                    filep = path;
                    filepSplit = regexp(filep,'\','split');
                    while size(filepSplit,2)>1&& ~isfile(filep)
                        filepSplit = regexp(filep,'\','split');
                        filep =join(filep(2:end),'\');
                    end
                    path = filep;
                end
                if ~isfile(path)
                    
                    
                    
                    path = [pwd,'\',fileName,ext];
                    if isfile(path)
                        info = getRecodringDetails([fileName,'.',ext],pwd,fileName);
                        S.info = info;
                    end
                end
                if ~flags %Dont load the vector, only the roi table
                    S.vec = [];
                    %                     info = getRecodringDetails([fileName,ext],pwd,fileName);
                    %                     S.info = info;
                    this = audioClip(S);
                    return
                end
                
                try
                    switch lower(S.sourceType)
                        
                        case 'wav'
                            vec = audioread(path) ;
                        case 'mat'
                            readfile = matfile(path);
                            classesInVar = whos(readfile);
                            indDouble = find(strcmpi('double',{classesInVar.class}),1,'first');
                            vec  = eval(['readfile.',classesInVar(indDouble).name]);
                        otherwise
                            error('Expected sourceType to be .mat or .wav but instead it was %s',S.sourceType)
                    end
                    S.vec = vec;
                    this         = audioClip(S);
                    
                    
                catch
                    
                    this = audioClip.loading_options();
                    
                end
            end
        end
    end
end

