function [net,info,inputs,outputs]=TrainMatchingNNAsMetric(Atlas,iterations,useParallel,trainingMode)
[P,G] = size(Atlas);        % P=number of atlas positions, G=number of atlas genes
inputs = zeros(2*G,P^2);    % each input is vector of position i and position j
outputs = zeros(P^2,1);     % output 1 if positions match, 0 if not

% set up parameters
trainFrac = 0.9;            % how much of input data to use as training data
noiseLevel = 0.1;          % level of noise to add to training data
noiseProb = .25;            % probability noise will be added during each 
                            % iteration of training
learningRate = 0.01;        % learning rate used in training (ADAM)
validationNumber = 20;     % how many times validation should be performed
validationPatience = Inf;    % stop if this many validation steps are
                            % worse than previous best

% default to no parallel
if nargin<3 || ~useParallel
    execEnvStr='auto';
else
    execEnvStr='parallel';
end

% default to trainingMode=1 (random split)
if nargin<4
    trainingMode = 1;
end

% build inputs/outputs matrices
currentIndex = 1;
for i=1:P
    for j=1:P
        inputs(:,currentIndex) = [Atlas(i,:)'; Atlas(j,:)'];
        outputs(currentIndex) = i==j;
        currentIndex = currentIndex+1;
    end
end

dataSize=size(inputs,1);        % should be 2*G
dataNumber=size(inputs,2);      % should be P^2

if trainingMode == 1
    % Training scheme 1: just split randomly
    % ----------
    trainNumber = ceil(trainFrac*dataNumber);
    testNumber = dataNumber-trainNumber;
    disp(['Out of ',num2str(dataNumber),' inputs, ',num2str(trainNumber), ...
        ' being used for training, ', num2str(testNumber),' for testing'])
    
    % split into train and test/validation set
    inputs = reshape(inputs, [1 1 dataSize dataNumber]);    % shape needed for trainNetwork
    testIndices = randperm(dataNumber,testNumber);  % random subset
    testInputs = inputs(:,:,:,testIndices);
    testOutputs = outputs(testIndices);
    trainInputs = inputs; trainInputs(:,:,:,testIndices)=[];
    trainOutputs = outputs; trainOutputs(testIndices)=[];

else
    % Training scheme 2: smaller test set with more 1's
    % ----------
    % split into train and test/validation set
    trainingMultiple = 99;
    % training set will include trainFrac% of inputs corresponding to outputs
    % of 1 along with this many times as many inputs with 0 output. For example,
    % if trainingMultiple=4, there will be 4x as many 0 outputs as 1 outputs,
    % so 1's will make up 1/(4+1)=1/5 of the training dataset
    exactMatches = find(outputs==1);                                        % how many 1's are in the dataset
    trainOnes = ceil(trainFrac*length(exactMatches));                       % how many 1's to train with
    trainIndices = exactMatches(randperm(length(exactMatches),trainOnes));  % where the 1's are
    nonMatches = find(outputs==0);                                          % how many 0's are in the dataset
    trainZeros = min(trainOnes*trainingMultiple,length(nonMatches));        % how many 0's to train with
    trainIndices = [trainIndices; nonMatches(randperm(length(nonMatches),trainZeros))]; % where the 0's are

    % display size of training dataset
    trainNumber = length(trainIndices);
    testNumber = dataNumber-trainNumber;
    frac = round(trainNumber/dataNumber*100,2);
    disp(['Out of ',num2str(dataNumber),' inputs, ',num2str(trainNumber), ...
        ' (' num2str(frac,2) '%) being used for training, ', ...
        num2str(testNumber),' (' num2str(100-frac,2) '%) for testing'])

    % do the splitting
    inputs = reshape(inputs, [1 1 dataSize dataNumber]);    % shape needed for trainNetwork
    trainInputs = inputs(:,:,:,trainIndices);
    trainOutputs = outputs(trainIndices);
    testInputs = inputs; testInputs(:,:,:,trainIndices)=[];
    testOutputs = outputs; testOutputs(trainIndices)=[];
end


% network structure: 2*G -> G -> G -> 1
layers = [
    imageInputLayer([1 1 dataSize],'Normalization','none')
    gaussianNoiseAtlasComparisonLayer(noiseLevel,0,1,noiseProb)   % custom noise filter
%     atlasSubtractionLayer("subtract")
    fullyConnectedLayer(G,"Name","fc_1")
    sigmoidLayer("sig_1")                   % custom activation function
    fullyConnectedLayer(G,"Name","fc_2")
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
    'Verbose',true,'Plots','training-progress');   

% train the network
[net,info]=trainNetwork(trainInputs,trainOutputs,layers,options);
end
