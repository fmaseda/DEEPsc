function [net,info,inputs,outputs]=TrainMatchingNNAsMetricUnlabelled(Atlas,iterations,useParallel)
[P,G] = size(Atlas);            % P=number of atlas positions, G=number of atlas genes
inputs = zeros(2*(G-1),G*P^2);  % each input is vector of position i and position j
                                % with one gene dropped out, total G*P^2
outputs = zeros(G*P^2,1);       % output 1 if positions match, 0 if not

% set up parameters
trainFrac = 0.9;            % how much of input data to use as training data
noiseLevel = 0.25;          % level of noise to add to training data
noiseProb = 1;              % probability noise will be added during each 
                            % iteration of training
learningRate = 0.01;        % learning rate used in training (ADAM)
validationNumber = 100;     % how many times validation should be performed
validationPatience = 50;    % stop if this many validation steps are
                            % worse than previous best

% default to no parallel
if nargin<3 || ~useParallel
    execEnvStr='auto';
else
    execEnvStr='parallel';
end

% build inputs/outputs matrices
currentIndex = 1;
for k=1:G       % dropout genes 1 by 1
    for i=1:P   % form "regular" inputs from Atlas without gene k
        for j=1:P
            inputs(:,currentIndex) = [Atlas(i,(1:G)~=k)'; Atlas(j,(1:G)~=k)'];
            outputs(currentIndex) = i==j;
            currentIndex = currentIndex+1;
        end
    end
end
dataSize=size(inputs,1);        % should be 2*(G-1)
dataNumber=size(inputs,2);      % should be G*P^2
trainNumber = ceil(trainFrac*dataNumber);
testNumber = dataNumber-trainNumber;
disp(['Out of ',num2str(dataNumber),' inputs, ',num2str(trainNumber), ...
    ' being used for training, ', num2str(testNumber),' for testing'])

% split into train and test/validation set
inputs = reshape(inputs, [1 1 dataSize dataNumber]);   % shape needed for trainNetwork
testIndices = randperm(dataNumber,testNumber);  % random subset
testInputs = inputs(:,:,:,testIndices);
testOutputs = outputs(testIndices);
trainInputs = inputs; trainInputs(:,:,:,testIndices)=[];
trainOutputs = outputs; trainOutputs(testIndices)=[];

% network structure: 2*(G-1) -> G-1 -> G-1 -> 1
layers = [
    imageInputLayer([1 1 dataSize],'Normalization','none')
    gaussianNoiseAtlasComparisonLayer(noiseLevel,0,1,noiseProb)   % custom noise filter
    fullyConnectedLayer(G-1,"Name","fc_1")
    sigmoidLayer("sig_1")                   % custom activation function
    fullyConnectedLayer(G-1,"Name","fc_2")
    sigmoidLayer("sig_2")
    fullyConnectedLayer(1,"Name","fc_3")
    sigmoidLayer("sig_out")
    ssePenaltyRegressionLayer("loss")];     % custom loss function

% training options
options = trainingOptions('adam', ...
    'InitialLearnRate',learningRate,...
    'Shuffle','every-epoch', ...            % shuffle training data each epoch
    'MiniBatchSize',trainNumber, ...
    'MaxEpochs',iterations, ...
    'ValidationData',{testInputs,testOutputs}, ...
    'ValidationFrequency',ceil(iterations/validationNumber), ...
    'ValidationPatience',validationPatience, ...
    'ExecutionEnvironment',execEnvStr, ...  % whether or not to use parallel
    ... % don't output training updates to command window, only plots
    'Verbose',false,'Plots','training-progress');   

% train the network
[net,info]=trainNetwork(trainInputs,trainOutputs,layers,options);
end
