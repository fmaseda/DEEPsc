function [net,info,inputs,outputs]=TrainDEEPsc(Atlas,varargin)
%% Trains a DEEPsc network on the provided spatial reference atlas
%
% Inputs
%   Atlas:              PxG array of gene expression data, where P is the number
%                       of atlas positions and G is the number of genes in the
%                       atlas
%
% Optional inputs
%   iterations:         integer, number of iterations to train (default=5000)
%   trainFrac:          0-1, fraction of data to use for training (default=0.9)
%   noiseLevel:         0-1, level of Gaussian noise to add to training data
%                       (default=0.25)
%   noiseProb:          0-1, probability that Gaussian noise will be added
%                       during each iteration of training (default=1)
%   learningRate:       real, learning rate used in training (ADAM, default=.01)
%   validationNumber:   integer, how many times validation should be
%                       performed during training (default=20)
%   validationPatience: integer, stop if this many validation steps are
%                       worse than previous best (default=Inf)
%   useParallel:        boolean, whether or not to use parallel processing
%                       if available (default=false)
%   trainingMode:       1 = uses all data and splits randomly, grows
%                       quadratically with size of Atlas;
%                       2 (default) = uses only fraction of matching targets,
%                       trainingMultiple times as many non-matching targets.
%   trainingMultiple:   integer, multiple of match targets to include as
%                       non-match targets in training data if trainingMode=2
%                       (default=99).
%   doPCA:              boolean, whether or not to perform PCA on Atlas
%                       before training (default=true)
%   PCAdims:            integer, number of PCA components to keep (default=8)
%
%
% Outputs
%   net:                The trained DEEPsc network
%   info:               Training info (loss, RMSE per iteration)
%   inputs:             4D array of training data, format that DEEPsc
%                       network requires
%   outputs:            array of labels for each input
%
% ----------
% Example usage
%   net = TrainMatchingNNAsMetric(MyAtlas,'iterations',5000,'useParallel',true)
%% Parse input arguments
for k = 1:2:length(varargin)
    switch lower(varargin{k})
        case 'iterations'
            iterations = varargin{k+1};
        case 'trainfrac'
            trainFrac = varargin{k+1};
        case 'noiselevel'
            noiseLevel = varargin{k+1};
        case 'noiseprob'
            noiseProb = varargin{k+1};
        case 'learningrate'
            learningRate = varargin{k+1};
        case 'validationnumber'
            validationNumber = varargin{k+1};
        case 'validationpatience'
            validationPatience = varargin{k+1};
        case 'useparallel'
            useParallel = varargin{k+1};
        case 'trainingmode'
            trainingMode  = varargin{k+1};
        case 'trainingmultiple'
            trainingMultiple = varargin{k+1};
        case 'dopca'
            doPCA = varargin{k+1};
        case 'pcadims'
            PCAdims = varargin{k+1};
    end
end

%% set defaults
if ~exist('iterations','var') || isempty(iterations)
    iterations = 5000;
end
if ~exist('trainFrac','var') || isempty(trainFrac)
    trainFrac = 0.9;
end
if ~exist('noiseLevel','var') || isempty(noiseLevel)
    noiseLevel = 0.25;
end
if ~exist('noiseProb','var') || isempty(noiseProb)
    noiseProb = 1;
end
if ~exist('learningRate','var') || isempty(learningRate)
    learningRate = .01;
end
if ~exist('validationNumber','var') || isempty(validationNumber)
    validationNumber = 20;
end
if ~exist('validationPatience','var') || isempty(validationPatience)
    validationPatience = Inf;
end
if ~exist('useParallel','var') || isempty(useParallel)
    useParallel = false;
end
if ~exist('trainingMode','var') || isempty(trainingMode)
    trainingMode = 2;
end
if ~exist('trainingMultiple','var') || isempty(trainingMultiple)
    trainingMultiple = 99;
end
if ~exist('usePCA','var') || isempty(doPCA)
    doPCA = true;
end
if ~exist('PCAdims','var') || isempty(PCAdims)
    PCAdims = 8;
end

% set up string for using parallel computing
if ~useParallel
    execEnvStr='auto';
else
    execEnvStr='parallel';
end
%% run PCA on Atlas
if doPCA
    [~,Atlas]=pca(Atlas*1);
    Atlas=Atlas(:,1:PCAdims);   % keep only PCAdims principal components

    % normalize to [0,1]
    Atlas=Atlas-min(Atlas);
    Atlas=NormalizeRNAseq(Atlas,'linear');
end

[P,G] = size(Atlas);        % P=number of atlas positions, G=number of atlas genes
inputs = zeros(2*G,P^2);    % each input is vector of position i and position j
outputs = zeros(P^2,1);     % output 1 if positions match, 0 if not

%% build inputs/outputs matrices
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
    % split into train and test/validation set;
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

%% form network
% network structure: 2*G -> G -> G -> 1
layers = [
    imageInputLayer([1 1 dataSize],'Normalization','none')
    gaussianNoiseAtlasComparisonLayer(noiseLevel,0,1,noiseProb)   % custom noise filter
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

%% train the network
[net,info]=trainNetwork(trainInputs,trainOutputs,layers,options);
end