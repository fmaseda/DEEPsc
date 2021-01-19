function [nets,foldIndices,info,Atlases]=TrainDEEPscPredRep(Atlas,varargin)
%% Trains a group of DEEPsc networks to be used for determining the predictive
%  reproducibility on a given reference atlas
%
% Inputs
%   Atlas:              PxG array of gene expression data, where P is the number
%                       of atlas positions and G is the number of genes in the
%                       atlas
%
% Optional inputs
%   numFolds:           integer, number of folds in k-fold cross validation
%                       (default=G, number of genes in Atlas)
%   foldIndices:        cell array, length=numFolds, each element containing
%                       indices of genes to predict in each of the folds; if
%                       not specified, fold indices are chosen randomly
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
%   showPlot:           boolean, whether or not to display plot of training
%                       progress (default=true)
%   Verbose:            boolean, whether or not to display training
%                       progress in command window (default=true)
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
%   doUMAP:             boolean, wehter or not to perform UMAP on Atlas
%                       before training (default=false); if doPCA=true, PCA
%                       is performed first, then UMAP is performed on the
%                       PCA space
%   UMAPdims:           integer, number of UMAP components to keep (default=8)
%   UMAPneighbors:      integer, number of neighbors to consider in UMAP
%                       algorithm (default=30)
%
%
% Outputs
%   nets:               Cell array of the trained DEEPsc networks
%   foldIndices:        Cell array of length numFolds, each element containing
%                       indices of genes to be reconstructed with each network
%   info:               Cell array of training info (loss, RMSE per iteration)
%   Atlases:            If dim reduction was done, this contains the
%                       reduced atlases for each network
%
% ----------
% Example usage
%   [nets,ind] = TrainDEEPscPredRep(MyAtlas,'iterations',5000,'useParallel',true,'numFolds',5)
%% Find numFolds and foldIndices arguments if there, remove them, and set default
for k = 1:2:length(varargin)
    switch lower(varargin{k})
        case 'numfolds'
            numFolds = varargin{k+1};
            numFoldsIndex=k;
        case 'foldindices'
            foldIndices = varargin{k+1};
            foldIndicesIndex=k;
    end
end
if exist('numFoldsIndex','var')
    varargin(numFoldsIndex:numFoldsIndex+1)=[];
end
if exist('foldIndicesIndex','var')
    % correct index if numFolds was supplied before foldIndices
    if exist('numFoldsIndex','var') && numFoldsIndex<foldIndicesIndex
        foldIndicesIndex = foldIndicesIndex - 2;
    end
    varargin(foldIndicesIndex:foldIndicesIndex+1)=[];
end
[~,G] = size(Atlas);	% G=number of atlas genes
if ~exist('numFolds','var') || isempty(numFolds)
    numFolds = G;
end

%% If foldIndices not supplied, randomly split atlas into folds for cross validation
if ~exist('foldIndices','var')
    allInd=randperm(G);
    foldSize=ceil(G/numFolds);
    foldIndices={};
    while length(allInd)>=foldSize
        % each fold contains foldSize randomly selected genes
        foldIndices{end+1} = allInd(1:foldSize);
        allInd(1:foldSize) = [];
    end
    if isempty(allInd)
        if length(foldIndices)~=numFolds
            % if not able to nicely divide folds, decrease number
            warning(['Atlas size ' num2str(G) ' with ' num2str(numFolds) ...
                ' folds leaves at least one fold empty; decreasing to ' ...
                num2str(length(foldIndices)) ' folds.'])
            numFolds = length(foldIndices);
        end
    else
        % otherwise, leftover indices in last fold
        foldIndices{end+1} = allInd;
    end
end

% error if supplied indices don't work
if numFolds~=length(foldIndices)
    error('Supplied foldIndices should be a cell array of length numFolds')
end

%% train each network on all other folds

% initialize outputs
nets=cell(1,numFolds);
info=cell(1,numFolds);
Atlases=cell(1,numFolds);

for k=1:numFolds
    ind=1:G;
    ind(foldIndices{k})=[];
    thisAtlas=Atlas(:,ind);
    
    % train network without showing plot
    [nets{k},info{k},Atlases{k}]=TrainDEEPsc(thisAtlas,varargin{:});
end
