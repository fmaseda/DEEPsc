function Corr=RunMatchingAlgorithmsDropout(Atlas,SCD,method,normMat,NN,numIter,Patches,doPCA)

% set variables for use in all methods
C=size(SCD,1);      % numberCells
P=size(Atlas,1);    % numberPositions
G=size(Atlas,2);    % numberGenes

% default input values
if nargin<4
    normMat=eye(G);
end
if nargin<5
    NN=1;
end
if nargin<6
    numIter=1;
end
if nargin<7
    Patches={};     % will cause error if method uses Patches
end
if nargin<8
    doPCA=false;
end

Corr=zeros(C,P,G);

for k=1:G
    Corr(:,:,k)=RunMatchingAlgorithms(Atlas(:,1:G~=k),SCD(:,1:G~=k),...
        method,normMat,NN,numIter,Patches,doPCA);
    % normMat(1:G~=k,1:G~=k)
end
