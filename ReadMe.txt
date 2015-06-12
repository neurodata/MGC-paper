The draft is named as RDC.pdf in the draft folder.


TibsSimu2Dim is the Matlab code for the first numerical section in RDC.
TibsSimu2Noise is the Matlab code for the second numerical section in RDC.
canoncorr suppress the full-rank warning for the above two codes; but it is not really used to compare rankdCorr and dCorr.


They are for the independence test of

y=f(Cx)+eps,

where f includes 10 types of functions, x is d*n, C is 1*d transformation, eps is 1*n white noise.



The pdf and the figures include 
1) the 10 tests with respect to dimension from d=1:1000 without noise (a noisy version looks similar),
2) the 10 test with respect to Gaussian noise at d=500,

for comparison, the figures further include
3) the 10 test with respect to Gaussian noise at d=1 (which is similar to the Tibs comments)
4) the 10 test with respect to Gaussian noise at d=300,
which are not included in the pdf write-up yet.
