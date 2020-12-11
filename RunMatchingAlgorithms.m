function Corr=RunMatchingAlgorithms(method,Atlas,SCD,varargin)
%% Runs the specified mapping method to map SCD to Atlas, outputting correspondence
%  scores for each pair
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
%   Corr:               CxP array with elements giving correspondence score
%                       for each cell-position pair
%
% ----------
% Example usage
%   Corr = RunMatchingAlgorithms('deepsc',MyAtlas,SCD,'NN',DEEPscNet,'doPCA',true)

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
if ~exist('numIter','var') || isempty(numIter)
    numIter = 1;
end
if ~exist('doPCA','var') || isempty(doPCA)
    doPCA = false;
end
if ~exist('PCAdims','var') || isempty(PCAdims)
    PCAdims = 8;
end

%% setup
C=size(SCD,1);      % numberCells
P=size(Atlas,1);    % numberPositions
G=size(Atlas,2);    % numberGenes

if doPCA
    if G>PCAdims
        G=PCAdims;
    end
    [PCA_coefs,Atlas]=pca(Atlas*1);
    Atlas=Atlas(:,1:G);
    Atlas=Atlas-min(Atlas);
    Atlas=NormalizeRNAseq(Atlas,'linear');
    
    SCD=SCD*PCA_coefs;
    SCD=SCD(:,1:G);
    SCD=SCD-min(SCD);
    SCD=NormalizeRNAseq(SCD,'linear');
end


%% run method
Corr=zeros(C,P);

switch lower(method)
    case {'binary','achim'}
% ------------------------------------------------------------
% Binary method from Achim (https://www.nature.com/articles/nbt.3209)
% ------------------------------------------------------------
        NormSCD=NormalizeRNAseq(SCD,'linear');  % scale SCD to [0,1] by gene to match atlas 
        Spec = NormSCD./mean(NormSCD);          % specificity score
        Spec(isnan(Spec))=0;
        expCutoff = .2;                 % number of counts to make expression 'true'
        BinSCD = NormSCD>expCutoff;     % binarized single cell data
        BinAtlas = Atlas>expCutoff;     % binarized atlas

        % calculate correspondence score, columns=location, rows=cells
        T=Spec./(1+Spec);
        Corr=(BinSCD.*T)*BinAtlas';     % positive contribution if expressed in
                                        % both SCD and atlas for that voxel
                                
        Corr=Corr-(BinSCD.*T)*((1-BinAtlas)'); % negative contribution if 
                                        % expressed in SCD but not atlas
                                
                                        % could also penalize here for being
                                        % expressed in atlas but not SCD
                                        
                                        
    case {'satija','seurat'}
% ------------------------------------------------------------
% Satija method for Seurat, uses multivariate Gaussians
% ------------------------------------------------------------
		expCutoff = .2;             % number of counts to make SCD expression 'true'
		BinAtlas = Atlas>expCutoff; % binarized Atlas
        SCD=SCD';                   % code from Julia expects SCD to be GxC, not CxG
        SCD=NormalizeRNAseq(SCD,'logscale'); % also needs to be log-normalized

		% skip LASSO-based imputation to ameliorate noise

		% now calculate Gaussian Mixture models for each marker gene
		fracOn=mean(BinAtlas,1);        % fraction of Atlas positions expressing
										% each marker gene; serves as cutoff for
										% Gaussian mixtures below
		cutoffs=zeros(G,1);
        for i=1:G	% loop through marker genes
			cutoffs(i)=quantile(SCD(i,:),1-fracOn(i));
        end
		groupAssignments=SCD>=cutoffs;	% true=on group, false=off group
		% each of the following are [numberGenes,2]-arrays; first column=off
		% Gaussian stats, second=on stats for each marker gene
		Gweights=[1-sum(groupAssignments,2)/C ...
					sum(groupAssignments,2)/C];
                % +eps() in denominator below to take care of case when all are on or off
		Gmeans=[sum(SCD.*~groupAssignments,2)./(sum(~groupAssignments,2)+eps()) ...
					sum(SCD.*groupAssignments,2)./(sum(groupAssignments,2)+eps())];
		Gstds=zeros(G,2);
        for i=1:G
			% std() of off group then on
			Gstds(i,:)=[std(SCD(i,~groupAssignments(i,:))) ...
                                std(SCD(i,groupAssignments(i,:)))];
        end
        % fix for all on or off
        Gstds(isnan(Gstds)|Gstds==0)=eps();

		% refine group assignments by moving each cell to group with closest mean
        if numIter>0
            for j=1:numIter
				meandists=zeros(G,C,2);
				meandists(:,:,1)=abs(SCD-Gmeans(:,1));
				meandists(:,:,2)=abs(SCD-Gmeans(:,2));
				groupAssignments=meandists(:,:,1)>=meandists(:,:,2);	% T=on, F=off
				% recalculate weights, means, stds with new group assignments
				Gweights=[1-sum(groupAssignments,2)/C ...
							sum(groupAssignments,2)/C];
				Gmeans=[sum(SCD.*~groupAssignments,2)./(sum(~groupAssignments,2)+eps()) ...
							sum(SCD.*groupAssignments,2)./(sum(groupAssignments,2)+eps())];
				Gstds=zeros(G,2);
                for i=1:G
					% std() of off group then on
					Gstds(i,:)=[std(SCD(i,~groupAssignments(i,:))) ...
                                        std(SCD(i,groupAssignments(i,:)))];
                end
                Gstds(isnan(Gstds)|Gstds==0)=eps();
%{
% 				% test plotting
%                 for i=1:G   % for every gene
%                     if i~=2
%                         continue
%                     end
%                     figure(i)
%                     hold off
%                     histogram(SCD(i,:),20,'Normalization','pdf')
%                     % off distribution
%                     hold on; plot(0:.01:max(SCD(i,:)), ...
%                               Gweights(i,1)*normpdf((0:.01:max(SCD(i,:)))', ...
%                                                       Gmeans(i,1),Gstds(i,1)))
%                     % on distribution
%                     hold on; plot(0:.01:max(SCD(i,:)), ...
%                               Gweights(i,2)*normpdf((0:.01:max(SCD(i,:)))', ...
%                                                       Gmeans(i,2),Gstds(i,2)))
%                     % combined distribution
%                     hold on; plot(0:.01:max(SCD(i,:)), ...
%                               Gweights(i,1)*normpdf((0:.01:max(SCD(i,:)))', ...
%                                                       Gmeans(i,1),Gstds(i,1)) ...
%                              +Gweights(i,2)*normpdf((0:.01:max(SCD(i,:)))', ...
%                                                       Gmeans(i,2),Gstds(2,2)))
%                     legend('','on','off','sum','Location','northeastoutside')
%                     title(['Gene ' num2str(i) ', iteration ' num2str(j)])
%                     pause
%                 end
%}
            end
        end

		% evaluate Gaussian components for each cell in SCD
		mixProbs=zeros(G,C,2);	% [numberGenes,numberCells,2]-array;
							    % first page="off" probs, second="on"
		% probability from "off" Gaussian
		for i=1:G   % loop through all genes in Atlas
			mixProbs(i,:,1)=Gweights(i,1)*normpdf(SCD(i,:)',Gmeans(i,1),Gstds(i,1));
		end
		% probability from "on" Gaussian
		for i=1:G   %#ok loop through all genes in Atlas
			mixProbs(i,:,2)=Gweights(i,2)*normpdf(SCD(i,:)',Gmeans(i,2),Gstds(i,2));
        end
		% normalize over cells
		mixProbs=mixProbs./(sum(mixProbs,2)+eps());
        
        AllProbs=zeros(G,C,P);
		% [numberGenes,numberCells,numberPositions]-array
        for j=1:C
            for i=1:G
                AllProbs(i,j,:)=reshape(mixProbs(i,j,BinAtlas(:,i)+1),[1 P]);
            end
        end
        
		for i=1:C   % loop through all cells
			AllProbs(:,i,:)=log(AllProbs(:,i,:))-max(AllProbs(:,i,:))- ...
                log(sum(exp(AllProbs(:,i,:)-max(AllProbs(:,i,:)))));
		end
 		AllProbs(AllProbs<-9.2)=-9.2;

		% initial mapping
        for i=1:C   % loop through all cells
			Corr(i,:)=reshape(exp(sum(AllProbs(:,i,:))),[P 1]);
% 			Corr(i,:)=reshape(exp(-sqrt(-sum(AllProbs(:,i,:)))),[P 1]);
			Corr(i,:)=Corr(i,:)/sum(Corr(i,:));
        end

        
    case {'karaiskos','matthews','mcc','distmap'}
% ------------------------------------------------------------
% Karaiskos method, uses Matthews correlation coefficient
% ------------------------------------------------------------
		NormSCD=NormalizeRNAseq(SCD,"linear"); % scale SCD to [0,1] by gene to match atlas
		expCutoff = .2;         % number of counts to make SCD expression 'true'
		BinSCD = NormSCD>expCutoff;     % binarized single cell data
		BinAtlas = Atlas>expCutoff;     % binarized Atlas

		TP=BinAtlas*BinSCD';            % only nonzero if both Atlas and SCD are on
		TN=(1-BinAtlas)*(1-BinSCD');    % only nonzero if both are off
		FP=(1-BinAtlas)*BinSCD';        % nonzero if on in SCD, off in atlas
		FN=BinAtlas*(1-BinSCD');        % nonzero if off in SCD, on in atlas
		% Corr=Matthews correlation coefficient
		Corr=(TP.*TN-FP.*FN)./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));
		Corr(isnan(Corr))=0;
        Corr=Corr';     % make dimensions consistent with other methods

        
    case {'peng','spearman','rcc'}
% ------------------------------------------------------------
% Peng method, uses Spearman's rank correlation coefficient
% ------------------------------------------------------------
        % assign ranks for genes in each position in Atlas 
        [~,idx]=sort(Atlas,2);
        [~,AtlasRank]=sort(idx,2);
        % assign ranks for genes in each cell in SCD 
        [~,idx]=sort(SCD,2);
        [~,SCDRank]=sort(idx,2);
        % calculate Pearson correlation coefficient
		Corr=corr(AtlasRank',SCDRank')';

        
    case {'infnorm','inf-norm','inf'}
% ------------------------------------------------------------
% use infinity norm to compare (scaled) SCD to atlas
% ------------------------------------------------------------
        NormSCD=NormalizeRNAseq(SCD,'linear'); % scale SCD to [0,1] by gene to match atlas
        for i=1:C       % loop through all cells in SCD
           for j=1:P    % for each cell, compare to location in atlas
               Corr(i,j)=1-max(NormSCD(i,:)-Atlas(j,:));
           end
        end
        
        
    case {'2norm','2-norm','2','euclidean'}
% ------------------------------------------------------------
% use 2-norm/euclidean distance
% ------------------------------------------------------------
        NormSCD=NormalizeRNAseq(SCD,'linear'); % scale SCD to [0,1] by gene to match atlas 
        for i=1:C       % loop through all cells in SCD
           for j=1:P    % for each cell, compare to location in atlas
               Corr(i,j)=1-norm(NormSCD(i,:)-Atlas(j,:))/sqrt(G);
               % division by size makes output in [0,1];
           end
        end
        
        
    case {'percent','percentdiff','percentdifference','%diff','% diff','%'}
% ------------------------------------------------------------
% use percent difference
% ------------------------------------------------------------
        NormSCD=NormalizeRNAseq(SCD,'linear'); % scale SCD to [0,1] by gene to match atlas
        for i=1:C   % iterate through all cells in dataset
        % calculate percentage difference between each element of atlas and SCD
        % expression, sum them, then normalize to get a correspondence score
        % for each cell in dataset to each cell in atlas
            percentDiff=2*abs(Atlas-repmat(NormSCD(i,:),P,1))./ ...
                (Atlas+repmat(NormSCD(i,:),P,1));
            percentDiff(isnan(percentDiff))=1;
            Corr(i,:)=1-sum(percentDiff,2)/(2*G);
        end
        
        
    case {'convolution-2norm','convolution-2','convolution-euclidean'}
% ------------------------------------------------------------
% use 2-norm with spatial convolution of the atlas first
% ------------------------------------------------------------
        ConvAtlas=Atlas;
        for i=1:P
            thisNbrs=[];
            % identify atlas positions that share an edge with position i
            for j=1:P
                if j==i
                    continue    % obviously neighbor of itself
                end
                if length(intersect(Patches{i},Patches{j}))>1
                    thisNbrs=[thisNbrs j]; %#ok<*AGROW>
                end
            end
            if ~isempty(thisNbrs)   % if cell has neighbors, convolve
                ConvAtlas(i,:)=(sum(Atlas(thisNbrs,:))+ ...
                    Atlas(i,:)*length(thisNbrs))/(2*length(thisNbrs));
            end
        end
        NormSCD=NormalizeRNAseq(SCD,'linear'); % scale SCD to [0,1] by gene to match atlas
        for i=1:C       % loop through all cells in SCD
           for j=1:P    % for each cell, compare to location in atlas
               Corr(i,j)=1-norm(NormSCD(i,:)-ConvAtlas(j,:))/sqrt(G);
               % division by size makes output in [0,1];
           end
        end
        
        
    case {'convolution-infnorm','convolution-inf'}
% ------------------------------------------------------------
% use inf-norm with spatial convolution of the atlas first
% ------------------------------------------------------------
        ConvAtlas=Atlas;
        for i=1:P
            thisNbrs=[];
            % identify atlas positions that share an edge with position i
            for j=1:P
                if j==i
                    continue    % obviously neighbor of itself
                end
                if length(intersect(Patches{i},Patches{j}))>1
                    thisNbrs=[thisNbrs j];
                end
            end
            if ~isempty(thisNbrs)   % if cell has neighbors, convolve
                ConvAtlas(i,:)=(sum(Atlas(thisNbrs,:))+ ...
                    Atlas(i,:)*length(thisNbrs))/(2*length(thisNbrs));
            end
        end
        NormSCD=NormalizeRNAseq(SCD,'linear'); % scale SCD to [0,1] by gene to match atlas
        for i=1:C       % loop through all cells in SCD
           for j=1:P    % for each cell, compare to location in atlas
               Corr(i,j)=1-max(NormSCD(i,:)-ConvAtlas(j,:));
           end
        end
    
        
    case {'weighted','weightednorm','matrixnorm','norm','lmnn'}
% ------------------------------------------------------------
% use weighted norm; requires normMat input
% ------------------------------------------------------------
        NormSCD=NormalizeRNAseq(SCD,'linear'); % scale SCD to [0,1] by gene to match atlas
        for i=1:C       % loop through all cells in SCD
           for j=1:P    % for each cell, compare to location in atlas
			   diff=NormSCD(i,:)-Atlas(j,:);
               Corr(i,j)=1/(1+norm(normMat*diff'));
           end
        end
        
        
    case {'neuralnet','ann','nn','deepsc'} 
% ------------------------------------------------------------
% use neural net; requires NN input
% ------------------------------------------------------------
        NormSCD=NormalizeRNAseq(SCD,'linear'); % scale SCD to [0,1] by gene to match atlas
        
        for i=1:C
            inputs=zeros(2*G,P);
            for j=1:P
                inputs(:,j)=[Atlas(j,:)'; NormSCD(i,:)'];
            end
            Corr(i,:)=predict(NN,reshape(inputs,[1 1 2*G P]),'MiniBatchSize',P);
        end
        
        
    case {'siamese','siamesenn','siam','snn'} 
% ------------------------------------------------------------
% use Siamese neural net; requires NN input to be cell array
% ------------------------------------------------------------
        NormSCD=NormalizeRNAseq(SCD,'linear'); % scale SCD to [0,1] by gene to match atlas
        for i=1:C
            inputs=zeros(2*G,P);
            for j=1:P
                inputs(:,j)=[Atlas(j,:)'; NormSCD(i,:)'];
            end
            Corr(i,:)=extractdata(predictSiamese(NN{1},NN{2}, ...
                reshape(inputs(1:G,:),[1 1 G P]),reshape(inputs(G+1:end,:), ...
                [1 1 G P])));
        end
            
% ------------------------------------------------------------
% more cases here
% ------------------------------------------------------------

    % 
    
    %
    
    % 
    
    %
    
    %
    
end

%% post-processing

if ~any(strcmpi({'neuralnet','ann','nn','siamese','siamesenn','siam','snn'},method))
    Corr=Corr-min(Corr,[],2);
    Corr=Corr./(max(Corr,[],2)+eps());  % divide each heatmap by maximum to make at
                                        % least one position equal one for every cell
end

end