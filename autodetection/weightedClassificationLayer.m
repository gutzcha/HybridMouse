classdef weightedClassificationLayer < nnet.layer.ClassificationLayer
    
    properties
        % (Optional) Layer properties.
        
        % Layer properties go here.
        ClassWeights
        
    end
    
    methods
        function layer = weightedClassificationLayer(classWeights, name, classNames)
            % layer = weightedClassificationLayer(classWeights) creates a
            % weighted cross entropy loss layer. classWeights is a row
            % vector of weights corresponding to the classes in the order
            % that they appear in the training data.
            %
            % layer = weightedClassificationLayer(classWeights, name)
            % additionally specifies the layer name.
            
            % Set class weights
            if ~isvector(classWeights)
               error('Invalid input, Class weights must be a vector') 
            end
%             if size(classWeights,2)~=1||size(classWeights,1)~=size(layer.Classes,1)
%                error('Expected calss weights to be of size %dx%d but they were of size %dx%d instead',...
%                    size(layer.Classes,1),size(layer.Classes,2),size(classWeights,1),size(classWeights,2)) 
%             end
            
            layer.ClassWeights = classWeights;
            
            % Set layer name and classes
            if nargin > 1
                layer.Name = name;
            end
            
            if nargin >2
                %Make sure that the number of classes and the number of
                %weights are the same
                assert (numel(classWeights)==numel(classNames))
                layer.Classes = classNames;
            end
            
            % Set layer description
            layer.Description = 'Weighted cross entropy';
            
        end
        
        function loss = forwardLoss(layer, Y, T)
            % loss = forwardLoss(layer, Y, T) returns the weighted cross
            % entropy loss between the predictions Y and the training
            % targets T.
            % Find observation and sequence dimensions of Y
            [~, N, S] = size(Y);
            
            % Reshape ClassWeights to KxNxS
            W = repmat(layer.ClassWeights(:), 1, N, S);
            
            % Compute the loss
            
            loss = -sum(sum(sum( W.*T.*log(nnet.internal.cnn.util.boundAwayFromZero(Y)))))/(N*S);
        end
        
        
        function dLdY = backwardLoss(layer, Y, T)
            % dLdY = backwardLoss(layer, Y, T) returns the derivatives of
            % the weighted cross entropy loss with respect to the
            % predictions Y.
            % Find observation and sequence dimensions of Y
            [~, N, S] = size(Y);
            
            % Reshape ClassWeights to KxNxS
            W = repmat(layer.ClassWeights(:), 1, N, S);
            
            % Compute the derivative
            dLdY = -(W.*T./nnet.internal.cnn.util.boundAwayFromZero(Y))/N;
        end
    end
end
