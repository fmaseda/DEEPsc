function [predrep,precision,WeightedCorr,predrepArray]=CorrErrorUnknownResult(Corr,Atlas,SCD)
% Corr is CxPxG, each page Corr from dropping out gene i
% Atlas is PxG
% SCD is CxG

% normalize SCD
SCD=NormalizeRNAseq(SCD,'linear');
% use Corr to predict value of dropped out gene for each cell based on
% other genes in the atlas
predError=zeros(size(Corr,[1 3])); % CxG array
for i=1:size(Corr,3)
    pred=((Corr(:,:,i).^8)./(sum((Corr(:,:,i).^8),2)+eps()))*Atlas(:,i); % predicted gene i value for each cell
    predError(:,i)=abs(pred-SCD(:,i));  % error in predicting gene i for each cell in SCD
end

weights=mean(predError);        % mean error across cells for each gene
weights=weights./(sum(weights)+eps());  % normalize so they add to 1
% create finalized Corr as weighted average (by prediction error) of each
% of the Corr pages
WeightedCorr=zeros(size(Corr,[1 2]));   % CxP array
for i=1:size(Corr,3)
    WeightedCorr = WeightedCorr+weights(i)*Corr(:,:,i);
end
% predictive reproducibility is mean error over genes then cells
predrepArray=mean(predError,2);
predrep=mean(predrepArray);
% calculate precision using weighted corr
precision=mean(abs((1-sum(WeightedCorr,2))./(size(WeightedCorr,2)-1)));
end