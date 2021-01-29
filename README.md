# DEEPsc

**DEEPsc** is a system-adaptive, deep learning-based method to impute spatial information onto a scRNA-seq dataset from a given spatial reference atlas. The following pipeline walks through how to train and implement a DEEPsc network.

---

## Quick start: Test on follicle data

Since the murine follicle system has the smallest reference atlas with *G=8* genes and thus takes the least amount of time to train, we have included a plug-and-play test for this system in the `/tests` folder. To run this test, simply navigate to the DEEPsc folder in MATLAB and run `TestFollicle` from the command window.

The following procedures can be followed for other systems or to better understand this test case:

---

## Importing the data

### Spatial reference atlases

A spatial reference atlas is stored as a *P x G* array, where *P* is the number of spatial positions in the atlas, and *G* is the number of genes for which there is known spatial expression. Each row contains the expression levels of each of the *G* genes.

The `/atlas` subfolder includes several `.mat` files which contain the matrices used for each of the three biological systems under study: the zebrafish embryo ([Satija et al., 2015](https://www.nature.com/articles/nbt.3192)), the Drosophila embryo ([Karaiskos et al., 2017](https://science.sciencemag.org/content/358/6360/194)), and the murine hair follicle ([Joost et al., 2016](https://www.cell.com/fulltext/S2405-4712(16)30265-4)). To load the data, simply navigate to the `/atlas` subfolder in the MATLAB "Current Folder" selector and double click to load each `.mat` file. Alternatively, run the following from the MATLAB Command Window:

```
load('/atlas/SYSTEM.mat')
```
where `SYSTEM` is replaced by either `Follicle`, `Drosophila`, or `Zebrafish`. For each system, several atlases are imported, including a binary and a continuous version, clearly labelled, e.g. `Atlas_ZebrafishBinary`.

### scRNA-seq datasets

A scRNA-seq dataset is stored as a *C x G* array, where *C* is the number of cells in the dataset, each row corresponding to a single cell. The preprocessed arrays for each of the three biological systems are included in `.mat` files in the `/scRNAseq` subfolder and can be loaded in the same way as the spatial reference atlases above. Each system has a continuous and a binarized version of scRNA-seq data, clearly labelled, e.g. `SCD_DrosophilaBinary`. The full, non-binarized scRNA-seq datasets (including genes for which no spatial information exists) for each system are stored in arrays such as `FullSCD_Drosophila`.

## Training a DEEPsc network

A DEEPsc network accepts a low-dimensional feature vector corresponding to a single position in the spatial reference atlas along with a corresponding feature vector of the gene expression of a single cell and returns a likelihood the input cell originated from the input position.

To train a network in MATLAB, run the following
```
net = TrainMatchingNNAsMetric(MyAtlas,'iterations',5000,'useParallel',true)
```
where `MyAtlas` is the (continuous) reference atlas you want to train on. Training options such as whether or not to perform PCA and how many prinicpal components to keep, how much noise to include in the training process, how many iterations to run, what kind of validation to use, and whether or not to use parallel processing if available, are all described in the documentation of `TrainMatchingNNAsMetric()`.

## Mapping cells with DEEPsc or other baseline methods

Once a DEEPsc network is trained, you can evaluate the network for various cell-pair combinations to determine the spatial origin of each cell. This software package also allows you to map cells with several other existing methods, including [Seurat](https://www.nature.com/articles/nbt.3192), [DistMap](https://science.sciencemag.org/content/358/6360/194), and several baseline comparison methods such as using a 2-norm or inf-norm as a metric, as well as a [Large margin nearest neighbor](https://en.wikipedia.org/wiki/Large_margin_nearest_neighbor) approach.

To map cells with a given method, call
```
Corr = RunMatchingAlgorithms(method,MyAtlas,SCD)
```
where `method` can be one of several strings (e.g. `'Seurat'`, `'DistMap'`, `'DEEPsc'`) defining which method to use, `MyAtlas` is the relevant reference atlas, and `SCD` is the scRNA-seq data. To run DEEPsc on a given atlas, call the following
```
Corr = RunMatchingAlgorithms('deepsc',MyAtlas,SCD,'NN',DEEPscNet,'doPCA',true)
```
The output `Corr` is a *C x P* array where `Corr(i,j)` contains the likelihood that cell `i` in the scRNA-seq dataset originated from point `j` in the spatial reference atlas.

## Quantifying spatial mapping performance

We define a comprehensive measure of evaluating how well a given method maps scRNA-seq data to known spatial origins, called a performance score. This score contains three components that measure the accuracy, precision, and robustness of a method, respectively. To determine accuracy, we map simulated scRNA-seq data with each cell being an exact replica of a single position in the reference atlas so that we can treat it as having a "known" origin. To calculate these values for a given method, call
```
[acc,prec,rob,perf] = MeasureMatchingRobustness(method,MyAtlas)
```
where `method` and `MyAtlas` are defined in the previous section. The output `acc` is the accuracy, `prec` is the precision, `rob` is the robustness, and `perf` is the overall performance score. By default, this shows a graphical summary of the robustness calculation, which requires multiple calculations of accuracy and precision with various levels of random Gaussian noise added to the simulated scRNA-seq data. This can be disabled by setting `'showPlot'` to `false`. To determine the performance score for a DEEPsc network without showing a plot, call
```
[acc,prec,rob,perf] = MeasureMatchingRobustness('deepsc',MyAtlas,'NN',DEEPscNet,'doPCA',true,'showPlot',false)
```

## Displaying inferred spatial origins

This software package outputs from `RunMatchingAlgorithms()` a *C x P* array of the likelihood a given method assigns each cell as having originated from each position. Thus, we can visualize rows of this array as a heatmap over a standard diagram of the system under study. To do this, call
```
DisplayMatchResults(system,MyCorr,cellNum)
```
where `system` can be one of three strings: `'follicle'`, `'zebrafish'`, or `'drosophila'`, and `cellNum` is the cell for which a heatmap is desired. There are also several required name-value pairs required to be passed to this function, as well as several optional parameters, all described in the documentation.

## Determining predictive reproducibility

To quantify the performance of a given method on actual scRNA-seq data where the origin is unknown, we also define a quantity called the predictive reproducibility, determined from a k-fold cross validation scheme where we split the reference atlas into k different folds, each containing some number of genes, then map the cells using all but one fold. We then use the heatmap of correspondence scores for each method to reconstruct each gene in the dropped out fold. To do this, we use the `CalculatePredictiveReproducibility()` function. For example, to determine the predictive reproducibility of DistMap using 5-fold cross validation, call
```
[pred_rep, array, values] = CalculatePredictiveReproducibility('distmap',MyAtlas,SCD,'numFolds',5)
```
where the second output is a 1xC array of the cell-by-cell predictive reproducibility, and the third output is the reconstructed atlas obtained by k-fold CV.

Note that to determine the predictive reproducibility of a DEEPsc network, additional training is required since the mapping with only a subset of genes cannot be done with the original network. Indeed for k-fold CV, k separate networks must be trained. To train a DEEPsc network to determine predictive reproducibility, use
```
[nets, indices] = TrainDEEPscPredRep(MyAtlas,'numFolds',k)
```
then calculate predictive reproducibility by setting the `NNs` property to `nets` and `foldIndices` to `indices`.