classdef LabeledAudioDatastore < matlab.io.Datastore & ...
        matlab.io.datastore.MiniBatchable& ...
        matlab.io.datastore.Shuffleable
    
    properties
        DatastoreSignal
        DatastoreNoise
        SignalNoiseRasio = 0.5
        Labels
        MiniBatchSize
        SequenceLength = 6;
        SamplingRate = 2.5e5;
        SequenceDimension
        FileSize
        Overlap
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
        
        function ds = LabeledAudioDatastore(signal_ds,noise_ds)
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
            ds.MiniBatchSize = 4;
            ds.NumObservations = numObservations;
            ds.CurrentFileIndex = 1;
            ds.CurrentFileIndex_noise = 1;
            ds.FileSize = 2;
            ds.Overlap = 0;
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
        
        function data = preprocessData(ds,predictors,responses)
            % data = preprocessData(ds,predictors,responses) preprocesses
            % the data in predictors and responses and returns the table
            % data
            
            predictors = cell2mat(vertcat(predictors{:})');
            responses = cell2mat(vertcat(responses{:})');
            
        
            fs = ds.SamplingRate;
            file_size = ds.FileSize;
            overlap = ds.Overlap;
            image_width = 51;
            feature_flag = false;
            new_fs_in = [];
                        
            defined_filter = define_a_filter(11,8);
            [predictors,~,responses]=...
                prepare_file_for_model(predictors,defined_filter,file_size,overlap,image_width,feature_flag,fs,new_fs_in,responses);
                        
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
            file_list = ds.Datastore.Files;
            labels = ds.Labels;
            
            sub_file_list = file_list(inds);
            sub_labels = labels(inds);
            
            sub_ds = copy(ds);
            sub_ds.Datastore.Files = sub_file_list;
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

