classdef spliceLayer < nnet.layer.Layer

    properties
        inputdim
    end

    properties (Learnable)
        % (Optional) Layer learnable parameters.

        % No learnable propeties
    end
    
    methods
        function layer = spliceLayer(inputdim,name)
            
         %Set name
         if exist('name','var')             
          layer.Name = name;
         end
         layer.inputdim =inputdim;
         %Description
          layer.Description =...
              'Multiplication layer of two inputs';
          
          %Set number of inputs
          layer.NumInputs = 1;
        end
        
        function Z = predict(layer,X1)
            Z = X1(:,layer.inputdim,:,:);
        end

       
    end
end