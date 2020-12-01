function [robustness,performance,accuracyArray,precisionArray]=MeasureMatchingRobustness(Atlas,SCD,method,knownOrigin,normMat,NN,showplot,numIter,Patches,doPCA)

P=size(Atlas,1);    % numberPositions
G=size(Atlas,2);    % numberGenes

% default input values
if nargin<4
    knownOrigin=true;   % default to assuming labelled data
end
if nargin<5
    normMat=eye(G);     % if method='weighted', default to identity matrix
end
if nargin<6
    NN=1;               % will cause error if method='ann' and no NN supplied
end
if nargin<7
    showplot=true;      % show plot by default
end
if nargin<8
    numIter = 1;
end
if nargin<9
    Patches={};
end
if nargin<10
    doPCA=false;
end

% calculate accuracy and precision with no noise
if knownOrigin
    assert(all(Atlas==SCD,'all'),'Assuming known spatial origin but SCD does not match Atlas');
    Corr=RunMatchingAlgorithms(Atlas,SCD,method,normMat,NN,numIter,Patches,doPCA);
    [accuracy,precision]=CorrErrorKnownResult(Corr)
else
%     SCD=NormalizeRNAseq(SCD,'linear');  % normalize SCD to be in [0,1]
    Corr=RunMatchingAlgorithmsDropout(Atlas,SCD,method,normMat,NN,numIter,Patches,doPCA);
    [accuracy,precision]=CorrErrorUnknownResult(Corr,Atlas,SCD)
end

% configure parameters
numSteps=20;
numRunsEachStep=10;

% calculate accuracy and precision for various levels of gaussian noise in
% the atlas
accuracyArray=zeros(numRunsEachStep,numSteps+1);    % +1 because of no noise step
precisionArray=zeros(numRunsEachStep,numSteps+1);
accuracyArray(:,1)=accuracy;
precisionArray(:,1)=precision;
tic
parfor i=1:numSteps
    for j=1:numRunsEachStep
        fprintf('i=%d, j=%d\n',i,j)
        NoisySCD=AddNoise(SCD,i/numSteps,0,1);
        if knownOrigin
            NoisyCorr=RunMatchingAlgorithms(Atlas,NoisySCD,method,normMat,NN,numIter,Patches,doPCA);
            [accuracy,precision]=CorrErrorKnownResult(NoisyCorr);
        else
            NoisyCorr=RunMatchingAlgorithmsDropout(Atlas,NoisySCD,method,normMat,NN,numIter,Patches,doPCA);
            [accuracy,precision]=CorrErrorUnknownResult(NoisyCorr,Atlas,SCD);
        end
        accuracyArray(j,i+1)=accuracy;
        precisionArray(j,i+1)=precision;
    end
end
toc

% accuracy, precision, total error
a=mean(accuracyArray); p=mean(precisionArray);
da=std(accuracyArray); dp=std(precisionArray);
e=1/2*(a+p);
de=sqrt(da.^2+dp.^2);

if showplot
    % plot results
    % ----------
    x=((1:numSteps+1)-1)./numSteps;
    % accuracy error
    subplot(2,2,1)
    errorbar(x,a,da,'go','MarkerSize',3)
    title('Accuracy Error')
    ylim([0 1])
    set(gca,'fontsize',14)

    % precision error
    subplot(2,2,2)
    errorbar(x,p,dp,'ro','MarkerSize',3)
    title('Precision Error')
    ylim([0 1])
    set(gca,'fontsize',14)

    % performance
    subplot(2,2,3)
    errorbar(x,e,de,'b-')
    title('Total error/Robustness')
    xlabel('Level of Gaussian Noise')
    ylim([0 1])
    set(gca,'fontsize',14)

    if knownOrigin
        sgtitle(['Robustness of method ''' method ''' for labelled data'])
    else
        sgtitle(['Robustness of method ''' method ''' for unlabelled data'])
    end
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

if showplot
    % and plot robustness on performance score plot
    hold on
    plot(x,theory*ones(size(x)),'b:')
    xline(1-r,'b:','LineWidth',3);
    hold off
end

robustness = r^4;
performance=1-(a(1)+p(1)+robustness)/3;

if showplot
    % display analytics
    txt=sprintf(['Accuracy: %9.4f\n\nPrecision: %8.4f\n\nRobustness:' ...
        ' %7.4f\n\nPerformance: %6.4f'],[a(1) p(1) robustness performance]);
    pan=uipanel('Title','Method Analytics','FontSize',14,'Position',[.5 0 .5 .5]);
    uicontrol(pan,'Style','text','Units', 'norm','Position',[0.15 -0.15 1 1],...
        'HorizontalAlignment','Left','FontSize',16,'FontName','FixedWidth', ...
        'FontWeight','Bold','String',txt)
    round((a(1)+p(1))/2,4)
end
