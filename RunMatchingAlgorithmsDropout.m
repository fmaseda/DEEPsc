function Corr=RunMatchingAlgorithmsDropout(method,Atlas,SCD,varargin)
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
%   normMat:            if method='weighted', this is the matrix of the
%                       weighted norm
%   NN:                 if method='deepsc', this is the trained DEEPsc
%                       network from TrainMatchingNNAsMetric()
%   numIter:            if method='seurat', this is the number of
%                       iterations to refine GMM modelling (default=1)
%   Patches:            if method requires convolution of the reference
%                       atlas, this contains the patches defined by the
%                       atlas
%   doPCA:              boolean, whether or not to perform PCA on Atlas/SCD
%                       before mapping (default=false)
%   PCAdims:            integer, number of PCA components to keep (default=8)
%
%
% Outputs
%   Corr:               CxPxG array with elements giving correspondence score
%                       for each cell-position pair. Pages correspond to
%                       correspondence scores for each dropped out gene.
%
% ----------
% Example usage
%   Corr = RunMatchingAlgorithmsDropout('deepsc',MyAtlas,SCD,'NN',DEEPscNet_Dropout,'doPCA',true)

%% parse input arguments
for k = 1:2:length(varargin)
    switch lower(varargin{k})
        case 'normmat'
            normMat = varargin{k+1};
        case 'nn'
            NN = varargin{k+1};
        case 'numiter'
            numIter = varargin{k+1};
        case 'patches'
            Patches = varargin{k+1};
        case 'dopca'
            doPCA = varargin{k+1};
        case 'pcadims'
            PCAdims = varargin{k+1};
    end
end

%% set defaults
C=size(SCD,1);      % numberCells
P=size(Atlas,1);    % numberPositions
G=size(Atlas,2);    % numberGenes

if ~exist('normMat','var') || isempty(normMat)
    normMat=eye(G);     % if method='weighted', default to identity matrix
end
if ~exist('NN','var') || isempty(NN)
    NN=1;               % will cause error if method='deepsc' and no NN supplied
end
if ~exist('numIter','var') || isempty(numIter)
    numIter = 1;
end
if ~exist('Patches','var') || isempty(Patches)
    Patches={};         % will cause an error if method requires it
end
if ~exist('doPCA','var') || isempty(doPCA)
    doPCA = false;
end
if ~exist('PCAdims','var') || isempty(PCAdims)
    PCAdims = 8;
end


%% run for each dropped out gene
Corr=zeros(C,P,G);

for k=1:G
    Corr(:,:,k)=RunMatchingAlgorithms(method,Atlas(:,1:G~=k),SCD(:,1:G~=k),...
        'normMat',normMat,'NN',NN,'numIter',numIter,'Patches',Patches,...
        'doPCA',doPCA,'PCAdims',PCAdims);
    % normMat(1:G~=k,1:G~=k)
end
