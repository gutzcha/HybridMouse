classdef sequenceDatastore_2 < matlab.io.Datastore & ...
        matlab.io.datastore.MiniBatchable & ...
        matlab.io.datastore.Shuffleable
    
    properties
        Datastore
        Labels
        SequenceDimension
        MiniBatchSize
        numObservations
        augmentData = false;
        augmentScale = false;
        YAugmentation = [-10,0]; %Image augmentation
        XAugmentation = [-50,10]; %Image augmentation
        AdditionRange = [-2 1]; %Image scale shift
        ScalenRange = 2; %Image scale shift
        AddLines = false;
        outputType = 'vector';
        windowSize = 51;
        noiseAddition = struct("speckle",     struct("flag",false,"value",0.01),...
            "gaussian",struct("flag",false,"value",  0.001));
        
        
    end
    
    properties(SetAccess = protected)
        NumObservations
        % This property is inherited from Datastore
        CurrentFileIndex
        
        outputOptions = {'map','vector'};
    end
    
    properties(Dependent)
        CurrentFileName
    end
    
    
    methods
        
        function ds = sequenceDatastore_2(folder)
            % Construct a MySequenceDatastore object
            
            % Create a file datastore. The readSequence function is
            % defined following the class definition.
            fds = datastore(...
                folder,'Type',"file","ReadFcn",@readSequence);
            
            ds.Datastore = fds;
            numObservations = numel(fds.Files);
            % Read labels from folder names
            %             for i = 1:numObservations
            %                 file = fds.Files{i};
            %                 filepath = fileparts(file);
            %                 [~,label] = fileparts(filepath);
            %                 labels{i,1} = label;
            %             end
            %             ds.Labels = categorical(labels);
            ds.Labels = categorical(1:numObservations);
            
            
            % Determine sequence dimension. When you define the LSTM
            % network architecture, you can use this property to
            % specify the input size of the sequenceInputLayer.
            X = preview(fds);
            X = X{1,1}{:,:};
            ds.SequenceDimension = size(X,1);
            
            % Initialize datastore properties.
            ds.MiniBatchSize = 5;
            ds.NumObservations = numObservations;
            ds.CurrentFileIndex = 1;
        end
        
        function tf = hasdata(ds)
            % Return true if more data is available
            tf = ds.CurrentFileIndex + ds.MiniBatchSize - 1 ...
                <= ds.NumObservations;
        end
        
        function [data,info] = read(ds)
            % Read one mini-batch batch of data
            miniBatchSize = ds.MiniBatchSize;
            info = struct;
            
            for i = 1:miniBatchSize
                ret = read(ds.Datastore);
                predictors(i,1)=ret{1,1};   %#ok
                responses(i,1)=ret{1,2};    %#ok
                vectors(i,1) = ret{1,3};    %#ok
                ds.CurrentFileIndex = ds.CurrentFileIndex + 1;
            end
            
            
            switch ds.outputType
                case 'map'
                    if ds.augmentData
                        [predictors,input] = ds.augmentImages(predictors,responses);
                    else
                        input = responses;
                    end
                    
                case 'vector'
                    
                    
                    %Change target label to one hot
                    %                     vectors = categorical2oh(vectors);
                    if ds.augmentData
                        [predictors,input] = ds.augmentImages(predictors,vectors);
                    else
                        input = vectors;
                    end
                otherwise
                    error('Invalid outputType, input must be %s',ds.outputOptions);
            end
            %Change window length
            if ds.windowSize>1
                predictors = expandWinWidth(ds,predictors);
            end
            
            %Return a table
            data = table(predictors,input);
        end
        
        function [predictors,input] = augmentImages(ds,predictors,input)
            % data = preprocessData(ds,predictors,responses) preprocesses
            % the data in predictors and responses and returns the table
            % data
            
                        data = cell2mat(predictors');
                        data = squeeze(permute(data,[4,1,2,3]));
            
                        %Add random lines
                        if ds.AddLines
                            xNum = 3;
                            yNum = 5;
                            sz = size(data,3);
                            for im = 1:sz
                                data(:,:,im) = addLines(data(:,:,im),xNum, yNum);
                            end
                        end
            
            
            
            
            miniBatchSize = ds.MiniBatchSize;
            for i = 1:miniBatchSize
                if ds.augmentData
                    X1 = squeeze(predictors{i});
                    if numel(size(X1))==2
                        X1 = permute(X1,[1,3,2]);
                    end
                    %The scale of the audio wave will affect the scale of
                    %the imput, here, I will randomly add values to the
                    %featre map to make the model robust to the scale and
                    %normalization of the audio imput. In anycase, the
                    %audio input MUST be zero centered
                    
                    %Rescale image 0-1
%                     X1 = mat2gray(X1);
                    
                    %Add random lines
                    if ds.AddLines
                        xNum = 3;
                        yNum = 5;
                        X1 = addLines(X1,xNum, yNum);
                    end
                    
                    %Add random noise                    
                    if ds.noiseAddition.speckle.flag
                        X1 = imnoise(X1,'speckle',ds.noiseAddition.speckle.value);
                    end
                    if ds.noiseAddition.gaussian.flag
                        
                        X1 = imnoise(X1,'gaussian',0,ds.noiseAddition.gaussian.value);
                    end
                    
                    %Change image scale and add bias
                    if ds.augmentScale
                        additionRange = ds.AdditionRange;
                        value_to_scale = max(normrnd(1,0.5),0.5);
                        value_to_add = (additionRange(2)-additionRange(1)).*rand(1,1) + additionRange(1);
                        X1 = X1.*value_to_scale+value_to_add;
                    end
                    
                    
                    fillVal = mean(X1,'all');
                    
                    X1_Cell = {X1};
                    
                    %Shift image up/down
                    if strcmp(ds.outputType,'map')
                        % Shift images vertiraclly.
                        %Get the images
                        augmenter = imageDataAugmenter( ...
                            'RandYTranslation',ds.YAugmentation,...
                            'RandXTranslation',ds.XAugmentation,...
                            'FillValue',fillVal);
                        
                        X2 = input{i};
                        X2_Cell = {X2};
                        X1_Cell = [X1_Cell,X2_Cell]; %#ok
                        augI = augment(augmenter,X1_Cell);
                        
                        %Remove added dimention
                        X2_aug = augI(1,2);
                        X2_aug = squeeze(X2_aug{:,:});
                        
                        T = adaptthresh(X2_aug);
                        X2_aug = imbinarize(X2_aug,T);
                        X2_aug = categorical(X2_aug);
                        X2_aug = permute(X2_aug,[4,2,1,3]);
                        
                        input{i} = X2_aug;
                        
                    else
                        X1_Cell = squeeze(X1_Cell{:});
                        
                        augmenter = imageDataAugmenter( ...
                            'RandYTranslation',ds.YAugmentation,...
                            'FillValue',fillVal);
                        augI = augment(augmenter,{X1_Cell});
                    end
                    %Add removed dimention
                    
                    X1_aug = augI(1,1);
                    X1_aug = X1_aug{:,:};
                    X1_aug = permute(X1_aug,[1,3,4,2]);
                    predictors{i} = X1_aug;
                end
                
            end
            function img = addLines(img,xNum, yNum)
                maxP = max(max(img)); %Get max value
                minP = min(min(img)); %Get min value;
                [numR,numC] = size(img); %Dimentions of image
                
                %pick number of vertical and horizontal lines randomly
                xNum = randi(xNum);
                yNum = randi(yNum);
                
                xLocs = randi(numR,[1,xNum]); %number of horizontal lines and their location
                yLocs = randi(numC,[1,yNum]); %number of vertical lines and their location
                
%                 randVec = @(len,maxVal) ones(1,len).*normrnd(maxVal*0.1,0.01,[1,len]); %Create vector with normal distributed random numbers values
                randVec = @(len,maxVal) ones(1,len).*normrnd(maxVal*0.5,0.01,[1,len]); %Create vector with normal distributed random numbers values
                
                %Add values to the image
                for indX = 1:xNum
                    vec = randVec(numC,maxP);
                    
                    rangerand = 1:randi(numC);
                    img(xLocs(indX),rangerand) = vec(rangerand);
                end
                
                for indY = 1:yNum
                    
                    
                    img(:,yLocs(indY)) = randVec(numR,maxP);
                end
                
                %Make sure that the added values don't pass the max val
                img(img>maxP) = maxP;
                img(img<minP) = minP;
                
                
                
            end
        end
        function reset(ds)
            % Reset to the start of the data
            reset(ds.Datastore);
            ds.CurrentFileIndex = 1;
        end
        function dsNew = shuffle(ds)
            % dsNew = shuffle(ds) shuffles the files and the
            % corresponding labels in the datastore.
            
            % Create a copy of datastore
            dsNew = copy(ds);
            dsNew.Datastore = copy(ds.Datastore);
            fds = dsNew.Datastore;
            
            % Shuffle files and corresponding labels
            numObservations = dsNew.NumObservations;
            idx = randperm(numObservations);
            fds.Files = fds.Files(idx);
            dsNew.Labels = dsNew.Labels(idx);
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
        function ind =get_CurrentFileIndex(ds)
            ind = ds.CurrentFileIndex;
        end
        function name = get.CurrentFileName(ds)
            name = ds.Datastore.Files(ds.CurrentFileIndex);
        end
        function predictors = expandWinWidth(ds,predictors)
            %Convert from 1 bin window to 51 bin window
            data = cell2mat(predictors');
            data = squeeze(permute(data,[4,1,2,3]));
            labels = [];
            defined_filter = [];
            %                 fs_model = 1997/6;
            file_size = []; %Sec - same as input
            overlap = 0; %Sec
            %                 image_width = 1; %Bins
            %                 image_width = 51; %Bins
            image_width = ds.windowSize;
            
            feature_flag = true; %Logical
            
            
            
            [predictors,~,~,~,~,~]=...
                prepare_file_for_model4(data,labels,defined_filter,file_size,overlap,image_width,feature_flag);
            
        end
    end
    
    methods (Hidden = true)
        
        function frac = progress(ds)
            % Determine percentage of data read from datastore
            frac = (ds.CurrentFileIndex - 1) / ds.NumObservations;
        end
        
    end
    
end % end class definition


function data = readSequence(filename)
% data = readSequence(filename) reads the sequence X from the MAT-file
% filename
mat_file = load(filename);
features = mat_file.features;
if isfield(mat_file,'labels_map')
    responses = mat_file.labels_map;
else
    responses = 0;
end

if isfield(mat_file,'labels_vector')
    vectors = mat_file.labels_vector;
else
    vectors = 0;
end


data={features,responses,vectors};
end
%
% function newImg = changeImgIntesity(img,sacle,offset)
% scaleVal = sacle(1) + (sacle(2)-sacle(1)) .* rand(1,1);
% offsetVal = offset(1) + (offset(2)-offset(1)) .* rand(1,1);
% img = img*scaleVal+offsetVal;
% end

function oh = categorical2oh(categoricalVector)
sz = size(categoricalVector,1);
oh = cell(sz,1);

for i = 1:sz
    tempVec = categoricalVector{i,1};
    ohTemp = zeros(2,length(tempVec));
    ohTemp(1,tempVec=='0') = 1;
    ohTemp(2,tempVec=='1') = 1;
    oh(i,1) = {categorical(ohTemp)};
end
end
