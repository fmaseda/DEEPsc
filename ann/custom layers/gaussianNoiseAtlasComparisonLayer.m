classdef gaussianNoiseAtlasComparisonLayer < nnet.layer.Layer
    % gaussianNoiseAtlasComparisonLayer   Gaussian noise layer
    %   A Gaussian noise layer adds random Gaussian noise pnly to the
    %   comparison part (lower half) of the Atlas input and thresholds it
    %   to fall in the interval[lowerlim, upperlim]
    %
    %   To create a Gaussian noise layer, use 
    %   layer = gaussianNoiseAtlasComparisonLayer(sigma, lowerlim, upperlim, name)

    properties
        Sigma       % Standard deviation.
        LowerLim    % lower limit after applying noise
        UpperLim    % upper limit after applying noise
        Prob        % probability that noise will be added during a
                    % single training iteration
    end
    
    methods
        function layer = gaussianNoiseAtlasComparisonLayer(sigma, lowerlim, upperlim, prob, name)
            % layer = gaussianNoiseAtlasComparisonLayer(sigma,lowerlim,upperlim,prob,name)
            % creates a Gaussian noise layer and specifies the standard
            % deviation, upper and lower cutoffs for and layer name.
            
            % default name = "noise"
            if nargin<5
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
            
            % default to 100% probability of adding noise
            if nargin<4
                layer.Prob = 1;
            else
                layer.Prob = prob;
            end
            
            layer.Sigma = sigma;
            layer.Description = ...
                "Gaussian noise on comparison cell with standard deviation " ...
                + layer.Sigma + " and threshold interval [" + layer.LowerLim ...
                + "," + layer.UpperLim + "], firing " + num2str(round(100*layer.Prob)) ...
                + "% of the time.";
            layer.Type = "Gaussian Noise";
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
            noise = randn(size(X)) * rand*sigma;    % extra rand makes added noise not always so large
            % don't add noise to first cell, only second
            noise(:,:,1:size(noise,3)/2)=0;
            Z=X;
            % only add noise Prob% of the time
            if rand<layer.Prob
                Z = X + noise;
                Z = max(Z,layer.LowerLim);
                Z = min(Z,layer.UpperLim);
            end
            
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