classdef sigmoidLayer < nnet.layer.Layer
    % sigmoid activation layer.
    
    methods
        function layer = sigmoidLayer(name) 
            % layer = sigmoidLayer(name) creates a sigmoid layer
            % and specifies the layer name.

            % Set layer name.
            if nargin<1
                layer.Name = "sigmoid";
            else
                layer.Name = name;
            end
            % Set layer description.
            layer.Description = "Sigmoid activation layer";
        end
        
        function Z = predict(layer, X)
            % Z = predict(layer, X) forwards the input data X through the
            % layer and outputs the result Z.
            
            Z = sigmoid(X);
        end
    end
end