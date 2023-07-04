classdef LabeledAudioDatastore < matlab.io.Datastore & ...
        matlab.io.datastore.MiniBatchable& ...
        matlab.io.datastore.Shuffleable
    
    properties
        DatastoreSignal
        DatastoreNoise
        SignalNoiseRasio = 0.95;
        Labels
        MiniBatchSize
        SequenceLength = 6;
        SamplingRate = 2.5e5;
        SequenceDimension
        FileSize
        Overlap
        ImageWidth
        OverrideExtractor
        extreactioOptionsForPreparation = struct(...
            'defined_filter',[],...
            'file_size',[],...
            'overlap',[],...
            'image_width',[],...
            'feature_flag',[],...
            'fs',[],...
            'new_fs_in',[],...
            'overrideExtractor',[])
        gaussina_filter = struct(...
            'win_len',10,...
            'g_std',2)
            
        
 
        
            
    end
    
    properties(Dependent)
        
         
    end
    
    properties(SetAccess = protected)
        NumObservations
    end
    
    properties(Access = private)
        % This property is inherited from Datastore
        CurrentFileIndex
        CurrentFileIndex_noise
    end
    
    
    methods
        
        function ds = LabeledAudioDatastore(signal_ds,noise_ds, extractor_param)
            % Construct a LabeledAudioDatastore object
   
            ds.DatastoreSignal = signal_ds;
            ds.DatastoreNoise = noise_ds;
            numObservations = numel(signal_ds.Files);
            ds.Labels = categorical(1:numObservations);
            % Determine sequence dimension. When you define the LSTM
            % network architecture, you can use this property to
            % specify the input size of the sequenceInputLayer.
            X = preview(ds.DatastoreSignal);
            ds.SamplingRate = X.sampling_rate;
            ds.SequenceDimension = X.numel_total;
            
            % Initialize datastore properties.
            ds.NumObservations = numObservations;
            ds.CurrentFileIndex = 1;
            ds.CurrentFileIndex_noise = 1;
            
            if ~exist('extractor_param','var')
               extractor_param = []; 
            end
              
            if isfield(extractor_param,'MiniBatchSize')
                ds.MiniBatchSize = extractor_param.MiniBatchSize;
            else
                ds.MiniBatchSize = 4;
            end
            
            if isfield(extractor_param,'FileSize')
                ds.FileSize = extractor_param.FileSize;
            else
                ds.FileSize = 2;
%                 ds.FileSize = 0.5;
            end
            % TODO: check input validity
            if isfield(extractor_param,'Overlap')
                ds.Overlap = extractor_param.Overlap;
            else
%                 ds.Overlap = 0;
                ds.Overlap = 5;
            end
            
            if isfield(extractor_param,'ImageWidth')
                ds.ImageWidth = extractor_param.ImageWidth;
            else
%                ds.ImageWidth = 51;
               ds.ImageWidth = 9;
            end
            
            if isfield(extractor_param,'OverrideExtractor')
                ds.OverrideExtractor = extractor_param.OverrideExtractor;
            else
                ds.OverrideExtractor = [];
            end
            
            
            
        end
        function set.gaussina_filter(ds,params)
            gaussian_params = ds.gaussina_filter;
            if isfield('win_len',params)
                gaussian_params.win_len = params.win_len;
            end
            if isfield('g_std',params)
                gaussian_params.g_std = params.g_std;
            end
            ds.gaussina_filter = gaussian_params;
        end
        function y = smooth_truth_vec(ds,x)
            gaussina_filter_param = ds.gaussina_filter;
            win_len = gaussina_filter_param.win_len;
            g_std = gaussina_filter_param.g_std;
%             w = gausswin(win_len, g_std)/sum(gausswin(win_len, g_std));
            w = gausswin(win_len, g_std);
            n_col = size(x,2);
            y = zeros(size(x));
            for ic = 1:n_col
                y(:,ic) = conv(x(:,ic),w','same');      
            end
            % clip
            y(y>1) = 1;
            
        end
        
        function tf = hasdata(ds)
            % Return true if more data is available
            tf = ds.CurrentFileIndex + ds.MiniBatchSize - 1 ...
                <= ds.NumObservations;
        end
        
          function reset(ds,type)
            % Reset to the start of the data
            if ~exist('type','var')||isempty(type)
               type = 'signal' ;
            end
            switch type
                case 'signal'
                    dsSource = ds.DatastoreSignal;
                    ds.CurrentFileIndex = 1;
                case 'noise'
                    dsSource = ds.DatastoreNoise;
                    ds.CurrentFileIndex_noise = 1;
                case 'both'
                    reset(ds,'signal')
                    reset(ds,'noise')
            end
            reset(dsSource);            
          end
        
        function [data,info] = read(ds)
            % Read one mini-batch batch of data
            
            %Decide from where to take, from signal to noise
            miniBatchSize = ds.MiniBatchSize;
            p = ds.SignalNoiseRasio;
            x = rand(1, miniBatchSize);
            source = (x<p);
            
            info = struct;
            
            for i = 1:miniBatchSize
                if source(i)
                    dsSource = ds.DatastoreSignal;
                    ret =  read(dsSource);
                else
                    dsSource = ds.DatastoreNoise;
                    if ~hasdata(dsSource)
                        reset(ds,'noise')
                    end
                    ret =  read(dsSource);
    
                end
                predictors{i,1} = ret.vector; %#ok
                responses(i,1) = {ret.labels}; %#ok
                ds.CurrentFileIndex = ds.CurrentFileIndex + 1; 
       
            end
            
            data = preprocessData(ds,predictors,responses);
        end
        
        function data = preprocessData(ds,signal,responses)
            % data = preprocessData(ds,predictors,responses) preprocesses
            % the data in predictors and responses and returns the table
            % data
            
            signal = cell2mat(vertcat(signal{:})');
%             responses = cell2mat(vertcat(responses{:})');
            
            % multilabel data
            responses_temp = (cat(2,responses{:}));
            responses_temp = cat(3, responses_temp{:});
            [n_bins, n_labels, n_vecs] = size(responses_temp);
            responses_cat = categorical(repmat("bg",[n_bins,n_vecs]));
            labels =  {'Low','subjectUSV','bg'};
            for icat = 1:n_labels-1
                for ivec = 1:n_vecs
                    responses_cat(squeeze(responses_temp(:,icat,ivec))) = labels{icat};
                end
            end
            
            %smooth responses
          
%             responses = smooth_truth_vec(ds,responses);
        
            fs = ds.SamplingRate;
            file_size = ds.FileSize;
            overlap = ds.Overlap;
            image_width = ds.ImageWidth;
            feature_flag = false;
            new_fs_in = [];
            output_class = 'categorical';
            override = ds.OverrideExtractor;
            
            %testing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
            
            windowDuration = round(fs*1); % one second
            overlapDuration = windowDuration-3; % 0.9 of window
            
            defined_filter = audioFeatureExtractor(...
                SampleRate=fs, ...
                Window=hann(windowDuration,"periodic"), ...
                OverlapLength=overlapDuration,...
                spectralCentroid=true);
%             override = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         signal = reshape(signal,[],1);
%                         if ~isempty(responses)
%                         responses = reshape(responses,[],1);
%                         end
%             defined_filter = define_a_filter(11,8);
            [predictors,~,responses,~,signal]=...
                prepare_file_for_model(signal,defined_filter,file_size,overlap,image_width,feature_flag,fs,new_fs_in,responses,override,output_class);
                        
           
            % Set option for future extraction
                        
            ret.defined_filter = defined_filter;
            ret.file_size = file_size;
            ret.overlap = overlap;
            ret.image_width = image_width;
            ret.feature_flag = feature_flag;
            ret.fs = fs;
            ret.new_fs_in = new_fs_in;
            ret.overrideExtractor = ds.OverrideExtractor;
            ds.extreactioOptionsForPreparation = ret;     
            
            % Return data as a table.
            data = table(predictors,responses);
        end
        
        function dsNew = shuffle(ds)
            % dsNew = shuffle(ds) shuffles the files and the
            % corresponding labels in the datastore.
            dsNew = copy(ds);
 
            fdsSignal = dsNew.DatastoreSignal;
                       
            numObservations = numel(fdsSignal.Files);
            idx = randperm(numObservations);
            fdsSignal.Files = fdsSignal.Files(idx);

            fdsNoise = dsNew.DatastoreNoise;
            
            numObservations = numel(fdsNoise.Files);
            idx = randperm(numObservations);
            fdsNoise.Files = fdsNoise.Files(idx);
           
        end
        function sub_ds = subset(ds,inds)
            if numel(inds)>ds.NumObservations||max(inds)>ds.NumObservations
                error('Indecies are out of bound')
            end
            file_list = ds.DatastoreSignal.Files;
            labels = ds.Labels;
            
            sub_file_list = file_list(inds);
            sub_labels = labels(inds);
            
            sub_ds = copy(ds);
            sub_ds.DatastoreSignal.Files = sub_file_list;
            sub_ds.Labels = sub_labels;
            sub_ds.NumObservations = numel(inds);
        end
        
    end
    
    methods (Hidden = true)
        
        function frac = progress(ds)
            % Determine percentage of data read from datastore
            frac = (ds.CurrentFileIndex - 1) / ds.NumObservations;
        end
        
    end
    
    
end % end class definition
function defined_filter = define_a_filter(win_len_ms,overlap_ms,spectType,fs)
% fs = 250000;
numBands = 124;

if ~exist("fs",'var')||isempty(fs)
    fs = 2.5e5;
end

w = round(fs*win_len_ms/1000);
overlap = round(fs*overlap_ms/1000);


if ~exist("spectType",'var')||isempty(spectType)
    spectType = 'melSpectrum';
    nfft =w*2;
end

if ~ischar(spectType)
    if spectType ==1
        spectType = 'melSpectrum';
        nfft =w*2;
    elseif spectType==2
        spectType = 'linearSpectrum';
        nfft = w;
    else
        error('Invalig type')
    end
end
% win_len_ms = 10;
% overlap_ms  =  6;


defined_filter = audioFeatureExtractor(...
    'Window',hamming(w,"periodic"),...
    'OverlapLength',overlap,...
    'SampleRate', fs,...
    'FFTLength',nfft,...
    'SpectralDescriptorInput',spectType,...
    spectType,true);

if strcmp(spectType,'melSpectrum')
    setExtractorParams(defined_filter,spectType,'NumBands',numBands,'FrequencyRange',[1000,125000]);
end
end

