The draft is named as RDC.pdf in the draft folder.

TibsSimu2Dim is the Matlab code for the first numerical section in RDC.
TibsSimu2Noise is the Matlab code for the second numerical section in RDC.
canoncorr suppress the full-rank warning for the above two codes; but it is not
really used to compare rankdCorr and dCorr.

They are for the independence test of

y=f(Cx)+eps,

where f includes 10+ types of functions, x is d*n, C is 1*d transformation, eps
is 1\*n noise.

The pdf and the figures include 1) the 10 tests with respect to dimension from
d=1:1000 without noise (a noisy version looks similar), 2) the 10 test with
respect to Gaussian noise at d=500,

for comparison, the figures further include 3) the 10 test with respect to
Gaussian noise at d=1 (which is similar to the Tibs comments) 4) the 10 test
with respect to Gaussian noise at d=300, which are not included in the pdf
write-up yet.

Results
=======

Theoretical Results
-------------------

1.  rdcorr test stat = 0 iff F(XY) = F(X) F(Y)

2.  dcorr is O(N\^2), HHG is O(N\^2 log N), rdcorr is ?  

Simulated Results
-----------------

### 1D Simulations

-   for 1D simulations designed to highlight HHG, n=30 & 50:

>   HHG \> rdcorr \> mdcorr, dcorr

-   for the other 10 1D functions, for n=30 & 50

rdcorr \\approx mdcorr,dcorr \> HHG

### n-D Simulations, n=50, d=1:500

-   for the nD settings, no noise,

rdcorr \approx mdcorr \> HHG

-   for the nD settings, multivariate noise, for n=50, d=1:500

rdcorr \>= mdcorr \> dcorr,HHG

### for the settings in which only rdcorr makes sense,

-   graphs and covariates

**rdcorr works?**

-   shapes and disease

rdcorr works?

Real Data Results
-----------------

-   in real data setting of CxP we have

p-value(rdcorr(k\*) \<\< p-value(rdcorr) \\approx p-value(HHG) \\approx 0.3 \<
p-value(mdcorr/dcorr)

-   in real data settings of shapes and diseases
