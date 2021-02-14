function [predrep,predrepArray_SCD,predrepArray_Atlas,predictedVals_SCD,predictedVals_Atlas,Corr] = ...
                CalculatePredictiveReproducibility(method,Atlas,SCD,varargin)
%% Runs the specified mapping method to map SCD to Atlas using a k-fold CV scheme,
%  outputting correspondence scores for each pair for each dropped out
%  gene, then reconstructing the dropped out gene based on the mapping
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
% Outputs
%   predrep:            length 4 array in [0,1], four measures of predictive
%                       reproducibility, [R_scd_zero,R_scd_nonzero,
%                                         R_atlas_zero,R_atlas_nonzero]
%   predrepArray_SCD:   Cx2 array of predictive reproducibility measures
%                       cell-by-cell
%   predrepArray_Atlas: Px2 array of predictive reproducibility measures
%                       position-by-position
%   predictedVals_SCD:  CxG array of reconstructed gene expression in SCD
%                       from cross validation
%   predictedVals_Atlas:PxG array of reconstructed gene expression in Atlas
%                       from cross validation
%   Corr:               cell array of CxP Corr matrices from each fold
%
% ----------
% Example usage
%   predrep = CalculatePredictiveReproducibility('2norm',MyAtlas,SCD,'numFolds',5)
%   [predrep, array_scd, array_atlas] = ...
%       CalculatePredictiveReproducibility('deepsc',MyAtlas,'numFolds',5, ...
%           'NNs',DEEPscPredRep_Nets,'FoldIndices', DEEPscPredRep_Indices, ...
%           'doPCA',true)

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

predrepArray_SCD=zeros(C,2);    % predictive reproducibility for each cell (1=zero, 2=nonzero)
predictedVals_SCD=zeros(C,G);    % reconstructed gene expression in SCD

predrepArray_Atlas=zeros(P,2);  % predictive reproducibility for each position (1=zero, 2=nonzero)
predictedVals_Atlas=zeros(P,G);  % reconstructed gene expression in Atlas

Corr=cell(1,numFolds);      % correspondence matrices used for reconstruction
predrep=zeros(G,4);         % [R_scd_zero,R_scd_nonzero,R_atlas_zero,R_atlas_nonzero]

for k=1:numFolds
    ind=1:G;
    ind(foldIndices{k})=[];
    % generate Corr using all other folds
    Corr{k}=RunMatchingAlgorithms(method,Atlas(:,ind),SCD(:,ind),...
        'normMat',normMat{k},'NN',NNs{k},'numIter',numIter,'Patches',Patches,...
        'doPCA',doPCA,'PCAdims',PCAdims);
    
    % reconstruct entire atlas with Corr, but only save this fold's genes
    C = Corr{k};
    C = C./sum(C);      % normalize over cells
    tmp = C'*SCD;
    predictedVals_Atlas(:,foldIndices{k}) = tmp(:,foldIndices{k});
    
    % same as above but now recontruct SCD
    C = Corr{k};
    C = C./sum(C,2);    % normalize over positions
    tmp = C*Atlas;
    predictedVals_SCD(:,foldIndices{k}) = tmp(:,foldIndices{k});
end

abs_diff_SCD = abs(predictedVals_SCD - SCD);
SCD_zero_mask = (SCD == 0);
abs_diff_Atlas = abs(predictedVals_Atlas - Atlas);
Atlas_zero_mask = (Atlas == 0);

% calculate predrep for each gene in SCD and in Atlas separately, then average over genes
predrep(:,1) = sum(abs_diff_SCD .* SCD_zero_mask,1) ./ sum(SCD_zero_mask,1);
predrep(:,2) = sum(abs_diff_SCD .* (1-SCD_zero_mask),1) ./ sum(1-SCD_zero_mask,1);
predrep(:,3) = sum(abs_diff_Atlas .* Atlas_zero_mask,1) ./ sum(Atlas_zero_mask,1);
predrep(:,4) = sum(abs_diff_Atlas .* (1-Atlas_zero_mask),1) ./ sum(1-Atlas_zero_mask,1);
predrep(isnan(predrep)) = 0;
predrep = 1-mean(predrep);	% mean over genes

% predictive reproducibility for zero and nonzero cells in SCD
predrepArray_SCD(:,1) = sum(abs_diff_SCD .* SCD_zero_mask,2) ./ sum(SCD_zero_mask,2);
predrepArray_SCD(:,2) = sum(abs_diff_SCD .* (1-SCD_zero_mask),2) ./ sum(1-SCD_zero_mask,2);
predrepArray_SCD(isnan(predrepArray_SCD)) = 0;

% predictive reproducibility for zero and nonzero positions in Atlas
predrepArray_Atlas(:,1) = sum(abs_diff_Atlas .* Atlas_zero_mask,2) ./ sum(Atlas_zero_mask,2);
predrepArray_Atlas(:,2) = sum(abs_diff_Atlas .* (1-Atlas_zero_mask),2) ./ sum(1-Atlas_zero_mask,2);
predrepArray_Atlas(isnan(predrepArray_Atlas)) = 0;

end
