classdef atlasSubtractionLayer < nnet.layer.Layer
    % atlasSubtractionLayer   Subtraction layer
    %   A subtraction layer takes as input a single vector with two parts
    %   corresponding to an input cell and an atlas location, subtracts the
    %   two, and produces a vector of half the original size
    %
    %   To create a subtraction layer, use 
    %   layer = atlasSubtractionLayer(name)
    
    methods
        function layer = atlasSubtractionLayer(name)
            % layer = atlasSubtractionLayer(name)
            % creates a subtaction layer and specifies the layer name.
            
            % default name = "subtract"
            if nargin<1
                layer.Name="subtract";
            else
                layer.Name = name;
            end
            
            layer.Description = ...
                "Subtraction layer.";
            layer.Type = "Subtraction";
        end
        
        function Z = predict(layer, X)
            % Z = predict(layer, X) forwards the input data X through the
            % layer for prediction and outputs the result Z.
            
            % Form matrix for autodifferentiable subtraction.
            G = size(X,3)/2;
            N = size(X,4);
            M = [eye(G) -eye(G)];
            Z = abs(M*squeeze(X));  % squeeze(X)
            Z = Z.reshape([1 1 G N]);
        end
        
    end
end