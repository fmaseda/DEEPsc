classdef gaussianNoiseWithThresholdLayer < nnet.layer.Layer
    % gaussianNoisewithThresholdLayer   Gaussian noise layer
    %   A Gaussian noise layer adds random Gaussian noise to the input and
    %   thresholds it to fall in the interval [lowerlim, upperlim]
    %
    %   To create a Gaussian noise layer, use 
    %   layer = gaussianNoiseWithThresholdLayer(sigma, lowerlim, upperlim, name)

    properties
        Sigma       % Standard deviation.
        LowerLim    % lower limit after applying noise
        UpperLim    % upper limit after applying noise
    end
    
    methods
        function layer = gaussianNoiseWithThresholdLayer(sigma, lowerlim, upperlim, name)
            % layer = gaussianNoiseWithThresholdLayer(sigma,name) creates a
            % Gaussian noise layer and specifies the standard deviation and
            % layer name.
            
            % default name = "noise"
            if nargin<4
                layer.Name="noise";
            else
                layer.Name = name;
            end
            
            % default to no thresholds
            if nargin<2
                layer.LowerLim = -Inf;
                layer.UpperLim = Inf;
            else
                layer.LowerLim = lowerlim;
                layer.UpperLim = upperlim;
            end
            
            layer.Description = ...
                "Gaussian noise with standard deviation " + sigma ...
                + "and threshold interval [" + lowerlim + "," ...
                + upperlim + "]";
            layer.Type = "Gaussian Noise";
            layer.Sigma = sigma;
        end
        
        function Z = predict(layer, X)
            % Z = predict(layer, X) forwards the input data X through the
            % layer for prediction and outputs the result Z.
            
            % At prediction time, the output is equal to the input.
            Z = X;
        end
        
        function [Z, memory] = forward(layer, X)
            % Z = forward(layer, X) forwards the input data X through the
            % layer and outputs the result Z.
            
            % At training time, the layer adds Gaussian noise to the input
            % and thresholds it.
            sigma = layer.Sigma;
            noise = randn(size(X)) * sigma;
            Z = X + noise;
            Z = max(Z,layer.LowerLim);
            Z = min(Z,layer.UpperLim);
            
            memory = [];
        end
        
        function dLdX = backward(layer, X, Z, dLdZ, memory)
            % [dLdX, dLdAlpha] = backward(layer, X, Z, dLdZ, memory)
            % backward propagates the derivative of the loss function
            % through the layer.
            
            % Since the layer adds a random constant, the derivative dLdX
            % is equal to dLdZ.
            dLdX = dLdZ;
        end
    end
end