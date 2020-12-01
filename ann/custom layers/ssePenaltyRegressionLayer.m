classdef ssePenaltyRegressionLayer < nnet.layer.RegressionLayer
    % regression layer with sum of squares error loss, with a larger
    % penalty if the network is supposed to output 1 than if it is supposed
    % to output 0
    
    methods
        function layer = ssePenaltyRegressionLayer(name)
            % layer = ssePenaltyRegressionLayer(name) creates a sum of
            % squares error regression layer and specifies the layer name.
    
            % Set layer name.
            if nargin<1
                layer.Name = "ssePenalty";
            else
                layer.Name = name;
            end

            % Set layer description.
            layer.Description = 'Sum of squares regression error with penalty';
        end
        
        function loss = forwardLoss(layer, Y, T)
            % loss = forwardLoss(layer, Y, T) returns the SSE loss between
            % the predictions Y and the training targets T.

            % Calculate sum of squares.
            sumSquares = sum(((Y-T).^2)./(1.001-T));
    
            % Take mean over mini-batch.
            N = length(Y);
            loss = sum(sumSquares)/N;
        end
    end
end