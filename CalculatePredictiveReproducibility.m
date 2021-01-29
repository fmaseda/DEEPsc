function [predrep,predrepArray,predictedVals]=CalculatePredictiveReproducibility(method,Atlas,SCD,varargin)
%% Runs the specified mapping method to map SCD to Atlas using a LOOCV scheme,
%  outputting correspondence scores for each pair for each dropped out gene
%
% Inputs
%   method:             string
%                       -'binary','achim': (Achim, 2015)
%                       -'satija','seurat': (Satija, 2015)
%                       -'karaiskos','matthews','mcc','distmap': (Karaiskos, 2017)
%                       -'peng','spearman','rcc': (Peng, 2016)
%                       -'infnorm','inf-norm','inf': Baseline
%                       -'2norm','2-norm','2','euclidean': Baseline
%                       -'percent','percentdiff','percentdifference','%diff','%diff','%': Baseline
%                       -'convolution-2norm','convolution-2','convolution-euclidean':
%                           2-norm baseline but convolved with spatial atlas,
%                           requires Patches to be supplied
%                       -'convolution-infnorm','convolution-inf': inf-norm 
%                           baseline but convolved with spatial atlas,
%                           requires Patches to be supplied
%                       -'weighted','weightednorm','matrixnorm','norm','lmnn':
%                           baseline weighted norm, requires normMat to be
%                           supplied
%                       -'neuralnet','ann','nn','deepsc': DEEPsc method,
%                           requires NN to be supplied
%                       -'siamese','siamesenn','siam','snn': Siamese neural
%                           network method, requies NN to be of format
%                           {mainNN, FClayer}
%   Atlas:              PxG array of gene expression data, where P is the number
%                       of atlas positions and G is the number of genes in the
%                       atlas
%   SCD:                CxG array of expression data, where C is the number
%                       of single cells, and G is the number of genes in
%                       the spatial reference atlas (indices should match
%                       Atlas)                           
%                       
%
% Optional inputs
%   normMat:            if method='weighted', this is a cell array of the
%                       weighted norms for each of the folds
%   NNs:                if method='deepsc', a cell array of the trained DEEPsc
%                       networks from TrainMatchingNNAsMetric(). Number of
%                       networks must match number of folds in
%                       cross-validation
%   numIter:            if method='seurat', this is the number of
%                       iterations to refine GMM modelling (default=1)
%   Patches:            if method requires convolution of the reference
%                       atlas, this contains the patches defined by the
%                       atlas
%   doPCA:              boolean, whether or not to perform PCA on Atlas/SCD
%                       before mapping (default=false)
%   PCAdims:            integer, number of PCA components to keep (default=8)
%   numFolds:           integer, number of folds in k-fold cross validation
%                       (default=G, number of genes in Atlas)
%   foldIndices:        cell array, length=numFolds, each element containing
%                       indices of genes to predict in each of the folds; if
%                       not specified, fold indices are chosen randomly
%
%
% Outputs
%   predrep:            double in [0,1], measure of predictive
%                       reproducibility
%   predrepArray:       length C array of predictive reproducibility
%                       cell-by-cell
%   predictedVals:      CxG array of reconstructed gene expression from
%                       cross validation
%
% ----------
% Example usage
%   predrep = CalculatePredictiveReproducibility('2norm',MyAtlas,SCD,'numFolds',5)
%   [predrep,array,vals] = CalculatePredictiveReproducibility('deepsc',MyAtlas, ...
%                       'numFolds',5,'NNs',{DEEPscNNFold1,DEEPscNNFold2, ...
%                       DEEPscNNFold3,DEEPscNNFold4,DEEPscNNFold5}, ...
%                       'FoldIndices', {DEEPscInd1,DEEPscInd2,DEEPscInd3, ...
%                       DEEPscInd4,DEEPscInd5},'doPCA',true)

%% parse input arguments
for k = 1:2:length(varargin)
    switch lower(varargin{k})
        case 'normmat'
            normMat = varargin{k+1};
        case {'nn','nns'}
            NNs = varargin{k+1};
        case 'numiter'
            numIter = varargin{k+1};
        case 'patches'
            Patches = varargin{k+1};
        case 'dopca'
            doPCA = varargin{k+1};
        case 'pcadims'
            PCAdims = varargin{k+1};
        case 'numfolds'
            numFolds = varargin{k+1};
        case 'foldindices'
            foldIndices = varargin{k+1};
        otherwise
            error(['Unrecognized input parameter ' varargin{k}])
    end
end

%% set defaults
C=size(SCD,1);      % numberCells
P=size(Atlas,1);    % numberPositions
G=size(Atlas,2);    % numberGenes

if ~exist('numIter','var') || isempty(numIter)
    numIter = 1;
end
if ~exist('Patches','var') || isempty(Patches)
    Patches={};     % will cause an error if method requires it
end
if ~exist('doPCA','var') || isempty(doPCA)
    doPCA = false;
end
if ~exist('PCAdims','var') || isempty(PCAdims)
    PCAdims = 8;
end
if ~exist('numFolds','var') || isempty(numFolds)
    numFolds = G;
end
if ~exist('NNs','var') || isempty(NNs)
    NNs=cell(1,numFolds);       % will cause error if method='deepsc' and no NNs supplied
end
if ~exist('normMat','var') || isempty(normMat)
    normMat=cell(1,numFolds);   % will cause an error if method requires it
end

%% split atlas into folds for cross validation
if ~exist('foldIndices','var')
    % if indices for folds not explicitly supplied, randomly split
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

%% Calculate predictive reproducibility for each fold

SCD=NormalizeRNAseq(SCD,'linear');  % normalize SCD
predrepArray=zeros(C,1);    % predictive reproducibility for each cell
predictedVals=zeros(C,G);   % reconstructed gene expression

for k=1:numFolds
    ind=1:G;
    ind(foldIndices{k})=[];
    Corr=RunMatchingAlgorithms(method,Atlas(:,ind),SCD(:,ind),...
        'normMat',normMat{k},'NN',NNs{k},'numIter',numIter,'Patches',Patches,...
        'doPCA',doPCA,'PCAdims',PCAdims);
    % use Corr to predict value of dropped out genes for each cell based on
    % other genes in the atlas
    for i=1:length(foldIndices{k})
        predictedVals(:,foldIndices{k}(i))=((Corr.^8)./(sum((Corr.^8),2)+eps())) ...
            *Atlas(:,foldIndices{k}(i));        % predicted gene i value for each cell
        predrepArray=predrepArray+abs(predictedVals(:,foldIndices{k}(i)) ...
            -SCD(:,foldIndices{k}(i)));         % error in predicting gene i for each cell in SCD
    end
end
predrepArray=1-predrepArray/G;    % mean error per gene
predrep=mean(predrepArray);     % single error value
