# Multiscale Generalized Correlation

- [repo contents](#repo-contents)
- [dependencies](#dependencies)
- [matlab](#matlab)
- [R](#R)

For pseudocode for all algorithms, see Appendix of draft in `Draft`.


## Repo Contents:

- [**Code**](https://github.com/neurodata-papers/MGC/tree/master/Code): folder containing MATLAB & R code to reproduce all results in the manuscript
- [**Draft**](https://github.com/neurodata-papers/MGC/tree/master/Draft): Discovering Relationships and their Structures Across Disparate Data Modalities,
Cencheng Shen, Carey E. Priebe, Mauro Maggioni, Qing Wang, Joshua T. Vogelstein,
submitted.
- [**Figures**](https://github.com/neurodata-papers/MGC/tree/master/Figures):  all figures from the plotting code used in the draft
- [**Data**](https://github.com/neurodata-papers/MGC/tree/master/Data):  contains the processed raw data to reproduce all results in the draft, and existing results to readily generate the figures.





## MATLAB

### Dependencies

We have tested in MATLAB R2017a on Windows 10 (on PC with i7 6850k and 64 GB memory) and  OSX Sierra (on MacBook Pro with 3.3 GHz Intel Core i7 and 16 GB RAM).

### Installation
Add all folders and subfolders of this repo to the path, you can do so using the command `addpath(genpath('<path-to-repo'))`, replacing `<path-to-repo>` with the path to your local copy of the repo.

### Demo
Type  `run_demo`
for a simulation example that outputs many things including visualization of the dependency, a p-value (pMGC) of < 0.05, the test statistic, and highlighted optimal scales in the multiscale significance map; it takes < 10 seconds to run.

### Reproduction Instruction

To reproduce the results in the manuscript, once installed, type any of the following:

- `plot_all` reproduce all figures in the draft from pre-generated results
- `run_1d_sims` runs the 1-dimensional simulations (~20 min)
- `run_hd_sims` runs high-dimensional simulations (~60 min)
- `run_realData` runs real data experiments (~20 min)


Note that the default number of replicates in each experiment is set at 100, which is much smaller than the number used in the draft. This can be increased by the function argument at the cost of linearly increasing the running time.

### Test on Real Data
1. To run on any given data X and Y, compute the n times n Euclidean distance matrices C for X and D for Y respectively, then type `MGCPermutationTest(C,D)`. If the input are already two distance matrices, use them directly.
2. The output will be the p-value, test statistic, and optimal scales. See the respective Matlab and R code for the output format.


## R

### Dependencies


We have tested in R v3.4.1 on Windows 10 (on PC with i7 6850k and 64 GB memory) and  OSX Sierra (on MacBook Pro with 3.3 GHz Intel Core i7 and 16 GB RAM).

On OSX, gfortan must also be installed, to do so, follow instructions from https://gcc.gnu.org/wiki/GFortranBinariesMacOS.


### installation

```
install.packages('ecodist')
install.packages('HHG')
install.packages('energy')
install.packages('SDMTools')
setwd("<path-to-repo>/Code/R")
source('MGCPermutationTest.R')
source('run_realData.R')
```

Note that `<path-to-repo>` must be replaced with the path to your local copy of the repo.



### Demo

Once installed, type `test=run_realData()` which takes < 1 minute to run.  The output is the variable test that includes:
- `localCorr` matrix
- `optimalInd` vector
- `pLocalCorr` matrix
- `pMGC` p-value
- `statMGC` test statistic.

Note that despite of the same implementation, the R version is slightly slower than Matlab for running on the same data, due to the fact that Matlab is slightly more efficient in handling matrix computation.
