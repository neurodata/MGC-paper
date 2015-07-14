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

 

1A) for 1D simulations designed to highlight HHG, n=30 & 50:

**HHG \> rdcorr \> mdcorr, dcorr**

1B) for the other 10 1D functions, for n=30 & 50

**rdcorr \\approx mdcorr,dcorr \> HHG**

2A) for the nD settings, no noise, for n=50, d=1:500

**rdcorr \\approx mdcorr \> HHG**

2B) for the nD settings, multivariate noise, for n=50, d=1:500

**rdcorr \>= mdcorr \> dcorr,HHG**

3A) for the settings in which only rdcorr makes sense, eg, graphs and covariates

**rdcorr works?**

3B) for another setting in which only rdocrr makes sense, eg, shapes and
covariates (not sure we need this)

**rdcorr works?**

4A) in real data setting of CxP we have

**p-value(rdcorr(k\*) \<\< p-value(rdcorr) \\approx p-value(HHG) \\approx 0.3 \<
p-value(mdcorr/dcorr)**

4B) in real data settings of shapes and diseases
