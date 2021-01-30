function [Norms, foldIndices]=TrainLMNNPredRep(Atlas,varargin)
%% Determines LMNN norms on a given reference atlas to calculate
% predictive reproducibility
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
%   doPCA:              boolean, whether or not to perform PCA on Atlas
%                       before training (default=true)
%   PCAdims:            integer, number of PCA components to keep (default=8)
%
%
% Outputs
%   Norms:              Cell array of the trained LMNN norms
%   foldIndices:        Cell array of length numFolds, each element containing
%                       indices of genes to be reconstructed with each matrix
%
% ----------
% Example usage
%   [Norms,ind] = TrainLMNNPredRep(MyAtlas,'numFolds',5)

%% Parse input arguments
for k = 1:2:length(varargin)
    switch lower(varargin{k})
        case 'numfolds'
            numFolds = varargin{k+1};
        case 'foldindices'
            foldIndices = varargin{k+1};
        case 'dopca'
            doPCA = varargin{k+1};
        case 'pcadims'
            PCAdims = varargin{k+1};
    end
end

%% set defaults
if ~exist('doPCA','var') || isempty(doPCA)
    doPCA = false;
end
if ~exist('PCAdims','var') || isempty(PCAdims)
    PCAdims = 50;
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

%% determine each metric using all other folds
Norms=cell(1,numFolds);
for k=1:numFolds
    disp(['Training metric on fold ' num2str(k) ' of ' num2str(numFolds) '...'])
    ind=1:G;
    ind(foldIndices{k})=[];
    thisAtlas=Atlas(:,ind);
    % dimensionality reduction for each fold
    if doPCA
        [~,thisAtlas]=pca(thisAtlas*1);
        thisAtlas=thisAtlas(:,1:PCAdims);
        % normalize to [0,1]
        thisAtlas=thisAtlas-min(thisAtlas);
        thisAtlas=NormalizeRNAseq(thisAtlas,'linear');
    end
    Norms{k}=lmnn2(thisAtlas);
end