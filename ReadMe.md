# Multiscale Generalized Correlation

- [repo contents](#repo-contents)
- [installation guide](#installation-guide)
- [reproduction instructions](#reproduction-instructions)

For pseudocode for all algorithms, see Appendix of draft in `Draft`.


## Repo Contents:

- [**Code**](https://github.com/neurodata-papers/MGC/tree/master/Code): folder containing MATLAB & R code to reproduce all results in the manuscript
- [**Draft**](https://github.com/neurodata-papers/MGC/tree/master/Draft): Discovering Relationships and their Structures Across Disparate Data Modalities, 
Cencheng Shen, Carey E. Priebe, Mauro Maggioni, Qing Wang, Joshua T. Vogelstein, 
submitted.
- [**Figures**](https://github.com/neurodata-papers/MGC/tree/master/Figures):  all figures from the plotting code used in the draft
- [**Data**](https://github.com/neurodata-papers/MGC/tree/master/Data):  contains the processed raw data to reproduce all results in the draft, and existing results to readily generate the figures.


## Installation guide:

### Dependencies

Either MATLAB or R, we have tested in MATLAB R2017a & R-3.4.1 on Windows 10. The local PC is equipped with i7 6850k and 64gb memory. 


### Test on Real Data
1. To run on any given data X and Y, compute the n times n Euclidean distance matrices C for X and D for Y respectively, then type `MGCPermutationTest(C,D)`. If the input are already two distance matrices, use them directly.
2. The output will be the p-value, test statistic, and optimal scales. See the respective Matlab and R code for the output format.


## Reproduction Instruction

### MATLAB

Add all folders and subfolders of MGC to the path. To reproduce all figures in the draft from pre-generated results, type
- `plot_all;` 

To re-generate results from scratch, type
- `run_1d_sims;`
- `run_hd_sims;`
- `run_realData;` 
which re-runs the 1-dimensional simulations, high-dimensional simulations, and real data experiments used in the draft. The running time is 20, 60, 20 minutes for each line above on the local PC.

Note that the default number of replicates in each experiment is set at 100, which is much smaller than the number used in the draft. This can be increased by the function argument at the cost of linearly increasing the running time.

### R

Set the working path to '/Code/R', and type 
- `test=run_realData`; 
to give the demo of MGC running on real data and output the results. It runs in 1 minute.



