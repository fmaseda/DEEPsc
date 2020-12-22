function varargout=MeasureMatchingRobustness(method,Atlas,varargin)
%% Determines performance score of method on Atlas
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
%
% Optional inputs
%   numSteps:           integer, number of steps between 0 and 1 to add
%                       noise at (default=20)
%   numRunsEachStep:    integer, number of runs to perform for each value
%                       of random noise added (default=10)
%   doPCA:              boolean, whether or not to perform PCA on Atlas/SCD
%                       before mapping (default=false)
%   PCAdims:            integer, number of PCA components to keep (default=8)
%   showPlot:           boolean, whether or not to show the statisticts
%                       in a new figure (default=true)
%   normMat:            if method='weighted', this is the matrix of the
%                       weighted norm
%   NN:                 if method='deepsc', this is the trained DEEPsc
%                       network from TrainMatchingNNAsMetric()
%   numIter:            if method='seurat', this is the number of
%                       iterations to refine GMM modelling (default=1)
%   Patches:            if method requires convolution of the reference
%                       atlas, this contains the patches defined by the
%                       atlas
%   useParallel:        boolean, whether or not to do calculations in
%                       parallel (default=true)
%
% Outputs
%   If nargout=1,       robustness
%   If nargout=2,       [robustness, performance]
%   If nargout=4,       [accuracy, precision, robustness, performance]
%   If nargout=6,       [accuracy, precision, robustness, performance, accArray, precArray]
%
% ----------
% Example usage
%   rob = MeasureMatchingRobustness('2norm',MyAtlas)
%   [rob,perf] = MeasureMatchingRobustness('deepsc',MyAtlas,'NN',DEEPscNet,'doPCA',true)
%   [acc,prec,rob,perf] = MeasureMatchingRobustness('lmnn',MyAtlas,'normMat',LMNN)
%   [acc,prec,rob,perf,accArray,precArray] = MeasureMatchingRobustness('inf',MyAtlas)

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
        case 'showplot'
            showPlot = varargin{k+1};
        case 'numsteps'
            numSteps = varargin{k+1};
        case 'numrunseachstep'
            numRunsEachStep = varargin{k+1};
        case 'useparallel'
            useParallel = varargin{k+1};
    end
end

%% default values
P=size(Atlas,1);    % numberPositions
G=size(Atlas,2);    % numberGenes

if ~exist('numSteps','var') || isempty(numSteps)
    numSteps=20;
end
if ~exist('numRunsEachStep','var') || isempty(numRunsEachStep)
    numRunsEachStep=10;
end
if ~exist('normMat','var') || isempty(normMat)
    normMat=eye(G);     % if method='weighted', default to identity matrix
end
if ~exist('NN','var') || isempty(NN)
    NN=1;               % will cause error if method='deepsc' and no NN supplied
end
if ~exist('showPlot','var') || isempty(showPlot)
    showPlot=true;      % show plot by default
end
if ~exist('numIter','var') || isempty(numIter)
    numIter = 1;
end
if ~exist('Patches','var') || isempty(Patches)
    Patches={};         % will cause an error if method requires it
end
if ~exist('doPCA','var') || isempty(doPCA)
    doPCA=false;
end
if ~exist('PCAdims','var') || isempty(PCAdims)
    PCAdims=8;
end
if ~exist('useParallel','var') || isempty(useParallel)
    useParallel=true;
end

%% begin calculation
% calculate accuracy and precision with no noise
Corr=RunMatchingAlgorithms(method,Atlas,Atlas,'normMat',normMat,'NN',NN, ...
        'numIter',numIter,'Patches',Patches,'doPCA',doPCA,'PCAdims',PCAdims);
[accuracy,precision]=CorrErrorKnownResult(Corr);

% calculate accuracy and precision for various levels of gaussian noise in
% the atlas
accuracyArray=NaN(numRunsEachStep,numSteps+1);    % +1 because of no noise step
precisionArray=NaN(numRunsEachStep,numSteps+1);
accuracyArray(:,1)=accuracy;
precisionArray(:,1)=precision;
tic
if useParallel
    parfor i=1:numSteps
        for j=1:numRunsEachStep
    %         fprintf('i=%d, j=%d\n',i,j)
            NoisyAtlas=AddNoise(Atlas,i/numSteps,0,1);
            NoisyCorr=RunMatchingAlgorithms(method,Atlas,NoisyAtlas,'normMat',normMat,'NN',NN, ...
                        'numIter',numIter,'Patches',Patches,'doPCA',doPCA,'PCAdims',PCAdims);
            [accuracy,precision]=CorrErrorKnownResult(NoisyCorr);
            accuracyArray(j,i+1)=accuracy;
            precisionArray(j,i+1)=precision;
        end
    end
else
    e=1/2*(accuracy+precision);
    theory=e+.1;
    i=1;
    while i<=numSteps && e<theory
        v1=accuracyArray(:,i+1);
        v2=precisionArray(:,i+1);
        for j=1:numRunsEachStep
    %         fprintf('i=%d, j=%d\n',i,j)
            NoisyAtlas=AddNoise(Atlas,i/numSteps,0,1);
            NoisyCorr=RunMatchingAlgorithms(method,Atlas,NoisyAtlas,'normMat',normMat,'NN',NN, ...
                        'numIter',numIter,'Patches',Patches,'doPCA',doPCA,'PCAdims',PCAdims);
            [accuracy,precision]=CorrErrorKnownResult(NoisyCorr);
            v1(j)=accuracy;
            v2(j)=precision;
        end
        accuracyArray(:,i+1)=v1;
        precisionArray(:,i+1)=v2;
        e=1/2*mean(v1+v2);
        i=i+1;
    end
end
toc % show time elapsed

% accuracy, precision, total error
a=mean(accuracyArray); p=mean(precisionArray);
da=std(accuracyArray); dp=std(precisionArray);
e=1/2*(a+p);
de=sqrt(da.^2+dp.^2);

if showPlot
    % plot results
    % ----------
    x=((1:numSteps+1)-1)./numSteps;
    % accuracy error
    subplot(2,2,1)
    errorbar(x,a,da,'go','MarkerSize',3)
    title('Accuracy Error')
    axis([0 1 0 1])
    set(gca,'fontsize',14)

    % precision error
    subplot(2,2,2)
    errorbar(x,p,dp,'ro','MarkerSize',3)
    title('Precision Error')
    axis([0 1 0 1])
    set(gca,'fontsize',14)

    % performance
    subplot(2,2,3)
    errorbar(x,e,de,'b-')
    title('Total error/Robustness')
    xlabel('Level of Gaussian Noise')
    axis([0 1 0 1])
    set(gca,'fontsize',14)

    sgtitle(['Robustness of method ''' method ''''])
end

% now calculate robustness
%theory=max(2/3,mean([2/3 performanceArray(1)]));
theory=e(1)+.1;
i=find(e>theory,1);
if isempty(i)
    r=0;
elseif i==1
    r=1;
else
    % 1-[amount of noise that causes performance score to drop below theory]
    r = theory+(i-2)*e(i)-e(i-1)*(i-1);
    r = 1 - r/(numSteps*(e(i)-e(i-1)));
end

if showPlot
    % and plot robustness on performance score plot
    hold on
    plot(x,theory*ones(size(x)),'b:')
    xline(1-r,'b:','LineWidth',3);
    hold off
end

robustness = r^4;
performance=1-(a(1)+p(1)+robustness)/3;
accuracy=a(1);
precision=p(1);

if showPlot
    % display analytics
    txt=sprintf(['Accuracy: %9.4f\n\nPrecision: %8.4f\n\nRobustness:' ...
        ' %7.4f\n\nPerformance: %6.4f'],[a(1) p(1) robustness performance]);
    pan=uipanel('Title','Method Analytics','FontSize',14,'Position',[.5 0 .5 .5]);
    uicontrol(pan,'Style','text','Units', 'norm','Position',[0.15 -0.15 1 1],...
        'HorizontalAlignment','Left','FontSize',16,'FontName','FixedWidth', ...
        'FontWeight','Bold','String',txt)
%     round((a(1)+p(1))/2,4)
end

if nargout==0||nargout==1
    varargout={robustness};
elseif nargout==2
    varargout={robustness,performance};
elseif nargout==4
    varargout={accuracy,precision,robustness,performance};
elseif nargout==6
    varargout={accuracy,precision,robustness,performance,accuracyArray,precisionArray};
end

end
