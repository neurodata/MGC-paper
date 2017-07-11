# Multiscale Generalized Correlation

This repo has:

- Code: 

folder containing MATLAB & R code to reproduce all results in the manuscript

- Draft: 

Discovering Relationships and their Structures Across Disparate Data Modalities, 
Cencheng Shen, Carey E. Priebe, Mauro Maggioni, Qing Wang, Joshua T. Vogelstein, 
submitted.

- Figures: 

all figures from the plotting code used in the draft

- Data: 

contains the processed raw data to reproduce all results in the draft, and existing results to readily generate the figures.



Installation guide:
In matlab, add all folders and subfolders of MGC to the path. 
To repeat the simulations and real data experiments, run:
>run_1d_sims;
>run_hd_sims;
>run_realData; 
To plot the figures, run:
>plot_all;
The running time on a standard i7 desktop takes around 1 day for 1D and HD simulations, and around 10 minutes for the real data. 


In R, source all codes in MGC/Code/R, and do
>run_realData
to give an example of MGC running on real data.
Typical running time: 1 minute