# Multiscale Generalized Correlation

## This repo has:

- **Code**: folder containing MATLAB & R code to reproduce all results in the manuscript
- **Draft**: Discovering Relationships and their Structures Across Disparate Data Modalities, 
Cencheng Shen, Carey E. Priebe, Mauro Maggioni, Qing Wang, Joshua T. Vogelstein, 
submitted.
- **Figures**: all figures from the plotting code used in the draft
- **Data**: contains the processed raw data to reproduce all results in the draft, and existing results to readily generate the figures.


## Installation guide:

### Dependencies

Either MATLAB or R, we have tested in MATLAB vXXX & R vXXX on OSX vXXX & ????

### Demo

1. To load demo data, type `load demo.mat`, which loads X and Y into the workspace
2. To run on data, in MATLAB, type `[a,b,c] = MGC(X,Y)`
3. The output will be a set of figures and the p-value, test statistic, optimal scales, XXX.




## Reproduction Instruction

### MATLAB

Add all folders and subfolders of MGC to the path. 
To repeat the simulations and real data experiments, run any of the following:
- `run_1d_sims;`
- `run_hd_sims;`
- `run_realData;` 
- `plot_all;` % to run all the plots

The running time on a standard i7 desktop takes around 1 day for 1D and HD simulations, and around 10 minutes for the real data. 

### R

All codes are in MGC/Code/R, and do `run_realData` to give an example of MGC running on real data.
Typical running time: 1 minute on a standard i7 desktop

