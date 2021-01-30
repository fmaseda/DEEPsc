%% load follicle data as example
load([pwd() '/atlases/Follicle.mat'])

%% train DEEPsc network
iterations = 20000;
DEEPscNet = TrainMatchingNNAsMetric(Atlas_FollicleContinuous,'iterations',iterations,'useParallel',true);

%% determine LMNN matrix for follicle atlas for baseline comparison
LMNN=lmnn2(Atlas_FollicleContinuous);

%% run DEEPsc along with various other mapping methods/baselines for comparison
Corr_Achim = RunMatchingAlgorithms('achim',Atlas_FollicleBinary,Atlas_FollicleContinuous);
Corr_Seurat = RunMatchingAlgorithms('seurat',Atlas_FollicleBinary,Atlas_FollicleContinuous);
Corr_Karaiskos = RunMatchingAlgorithms('distmap',Atlas_FollicleBinary,Atlas_FollicleContinuous);
Corr_Peng = RunMatchingAlgorithms('peng',Atlas_FollicleContinuous,Atlas_FollicleContinuous);
Corr_2norm = RunMatchingAlgorithms('2',Atlas_FollicleContinuous,Atlas_FollicleContinuous);
Corr_Infnorm = RunMatchingAlgorithms('inf',Atlas_FollicleContinuous,Atlas_FollicleContinuous);
Corr_Percent = RunMatchingAlgorithms('%',Atlas_FollicleContinuous,Atlas_FollicleContinuous);
Corr_LMNN = RunMatchingAlgorithms('lmnn',Atlas_FollicleContinuous,Atlas_FollicleContinuous, ...
    'normMat',LMNN);
Corr_DEEPsc = RunMatchingAlgorithms('deepsc',Atlas_FollicleContinuous,Atlas_FollicleContinuous, ...
    'NN',DEEPscNet,'doPCA',true);

%% compile performance statistics and display as table
stats=zeros(9,4);
[stats(1,1),stats(1,2),stats(1,3),stats(1,4)] = MeasureMatchingRobustness('achim',Atlas_FollicleBinary,'showPlot',false);
[stats(2,1),stats(2,2),stats(2,3),stats(2,4)] = MeasureMatchingRobustness('seurat',Atlas_FollicleBinary,'showPlot',false);
[stats(3,1),stats(3,2),stats(3,3),stats(3,4)] = MeasureMatchingRobustness('distmap',Atlas_FollicleBinary,'showPlot',false);
[stats(4,1),stats(4,2),stats(4,3),stats(4,4)] = MeasureMatchingRobustness('peng',Atlas_FollicleContinuous,'showPlot',false);
[stats(5,1),stats(5,2),stats(5,3),stats(5,4)] = MeasureMatchingRobustness('2',Atlas_FollicleContinuous,'showPlot',false);
[stats(6,1),stats(6,2),stats(6,3),stats(6,4)] = MeasureMatchingRobustness('inf',Atlas_FollicleContinuous,'showPlot',false);
[stats(7,1),stats(7,2),stats(7,3),stats(7,4)] = MeasureMatchingRobustness('%',Atlas_FollicleContinuous,'showPlot',false);
[stats(8,1),stats(8,2),stats(8,3),stats(8,4)] = MeasureMatchingRobustness('lmnn',Atlas_FollicleContinuous,...
                                                                            'normMat',LMNN,'showPlot',false);
[stats(9,1),stats(9,2),stats(9,3),stats(9,4)] = MeasureMatchingRobustness('deepsc',Atlas_FollicleContinuous,...
                                                                            'NN',DEEPscNet,'doPCA',true,'showPlot',false);
array2table(stats,...
    'RowNames',{'Achim','Seurat','DistMap','Peng','2-norm','Inf-norm','%diff','LMNN','DEEPsc'},...
    'VariableNames',{'accuracy','precision','robustness','performance'})
                                                                        
%% display mapping of random cell

i=randi(size(Atlas_FollicleContinuous,1));
DisplayMatchResults('follicle',Corr_Achim,i,'Patches',Patches_Follicle,...
    'Vertices',Vertices_Follicle,'title',['Cell ' num2str(i) ', Achim']);
DisplayMatchResults('follicle',Corr_Seurat,i,'Patches',Patches_Follicle,...
    'Vertices',Vertices_Follicle,'title',['Cell ' num2str(i) ', Seurat']);
DisplayMatchResults('follicle',Corr_Karaiskos,i,'Patches',Patches_Follicle,...
    'Vertices',Vertices_Follicle,'title',['Cell ' num2str(i) ', DistMap']);
DisplayMatchResults('follicle',Corr_Peng,i,'Patches',Patches_Follicle,...
    'Vertices',Vertices_Follicle,'title',['Cell ' num2str(i) ', Peng']);
DisplayMatchResults('follicle',Corr_2norm,i,'Patches',Patches_Follicle,...
    'Vertices',Vertices_Follicle,'title',['Cell ' num2str(i) ', 2-norm']);
DisplayMatchResults('follicle',Corr_Infnorm,i,'Patches',Patches_Follicle,...
    'Vertices',Vertices_Follicle,'title',['Cell ' num2str(i) ', Inf-norm']);
DisplayMatchResults('follicle',Corr_Percent,i,'Patches',Patches_Follicle,...
    'Vertices',Vertices_Follicle,'title',['Cell ' num2str(i) ', % Difference']);
DisplayMatchResults('follicle',Corr_LMNN,i,'Patches',Patches_Follicle,...
    'Vertices',Vertices_Follicle,'title',['Cell ' num2str(i) ', LMNN']);
DisplayMatchResults('follicle',Corr_DEEPsc,i,'Patches',Patches_Follicle,...
    'Vertices',Vertices_Follicle,'title',['Cell ' num2str(i) ', DEEPsc']);

%% Load scRNA-seq and compute predictive reproducibility
load([pwd() '/scRNAseq/Follicle.mat'])
% train new DEEPsc ensembles first
[DEEPscPredRepNets1,indices]=TrainDEEPscPredRep(Atlas_FollicleContinuous,'numFolds',4);
DEEPscPredRepNets2=TrainDEEPscPredRep(Atlas_FollicleContinuous,'numFolds',4,'foldIndices',indices);
DEEPscPredRepNets3=TrainDEEPscPredRep(Atlas_FollicleContinuous,'numFolds',4,'foldIndices',indices);

% train new LMNN metrics also
LMNNPredRep=TrainLMNNPredRep(Atlas_FollicleContinuous,'numFolds',4,'foldIndices',indices);

% compute pred. rep. for each method
predReps = zeros(9,1);
predReps(1) = CalculatePredictiveReproducibility('achim',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices);
predReps(2) = CalculatePredictiveReproducibility('seurat',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices);
predReps(3) = CalculatePredictiveReproducibility('distmap',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices);
predReps(4) = CalculatePredictiveReproducibility('peng',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices);
predReps(5) = CalculatePredictiveReproducibility('2',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices);
predReps(6) = CalculatePredictiveReproducibility('inf',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices);
predReps(7) = CalculatePredictiveReproducibility('%',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices);
predReps(8) = CalculatePredictiveReproducibility('lmnn',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices,'normMat',LMNNPredRep);

% for DEEPsc, average over three results
temp = CalculatePredictiveReproducibility('deepsc',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices,'NNs',DEEPscPredRepNets1);
temp = temp+CalculatePredictiveReproducibility('deepsc',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices,'NNs',DEEPscPredRepNets2);
temp = temp+CalculatePredictiveReproducibility('deepsc',Atlas_FollicleContinuous, ...
    SCD_Follicle_normalized,'numFolds',4,'foldIndices',indices,'NNs',DEEPscPredRepNets3);
predReps(9) = temp/3;

array2table(predReps,...
    'RowNames',{'Achim','Seurat','DistMap','Peng','2-norm','Inf-norm','%diff','LMNN','DEEPsc'},...
    'VariableNames',{'predictive reproducibility'})