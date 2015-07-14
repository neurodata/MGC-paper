title: Multimodal Two-Sample Testing authors: sigma(Cencheng Shen, cep, mauro?,
jovo)

Introduction
============

-   **O** - Information age, lots of data, lots of questions. Data are often
    high-dimensional or non-Euclidean. Often want to test whether one
    view/modality of an entity is independent of another.

-   **G** - There does not exist a two-sample testing methodology that can be
    applied to arbitrary metric data

-   **C** - 1D and nD methods exist. but nothing for non-euclidean data.

-   **A** - rankdcorr

-   **R** - our test is fast, easy to implement, has strong theoretical support,
    is state of the art, and works in settings for which no other known method
    does.

-   **F** - opens the door to address scientific questions that we previously
    intractable.

# Methods


## Rankdcorr


describe rankdcorr here

## Simulations


### Euclidean Simulations

They are for the independence test of

$$
y=f(Cx)+eps,
$$

where f includes 10+ types of functions, x is $d \times n$, C is $1\times d$
transformation, eps is 1/\*n noise.


### Graph vs Covariate Simulations

$n$ vertices per graph

$m$ subjects

$D$ is the maximum rank of the expected adjacency matrices

$X_{d} \in R^{n}$ of latent eigenvector

$\lambda_{i}$  $\in \triangle_{D}$ is the relative weights of the eigenvectors per subject i

$P_{i} \in [0,1]^{n x n}$ is the expected adjacency matrix for subject i

$A_{i}$  is the adjacency matrix for subject i, trial j

$Y_{i} \in R^p$ is vector of covariates

Some proposed constants:

$n=100$, $m=50$, $p=5$

X_d = e_d, that is, X_d is a zero vector with a single 1 in the d-th element

$Y \sim F_Y = Dir(1)$, so, uniform

$\lambda | y \sim F_{\lambda|Y} = Dir( Y )$, 



````
for all i in [n]
    $y_{i} \sim F_Y$
    $\lambda_{i} | y \sim$ F_{\lambda | y}
    $P_{i} = \sum_{d}^{D} \lambda_{id} \langle X_{d}, X_{d} \rangle$
    $A_{i} \sim$  Bern $(\langle X_{i}, X_{i} \rangle)$ 
end
````

Results
=======

Theoretical Results
-------------------

1.  rdcorr test stat = 0 iff F(XY) = F(X) F(Y)

2.  dcorr is O(N\^2), HHG is O(N\^2 log N), rdcorr is ?  

Simulated Results
-----------------

### 1D Simulations, n=30, 50

#### A) no noise

-   for 1D simulations designed to highlight HHG, n=30 & 50:

$$
HHG > rdcorr > mdcorr, dcorr
$$

-   for the other 10 1D functions, for n=30 & 50

$$
rdcorr \approx mdcorr,dcorr > HHG
$$

#### B) additive noise

?

#### C) multiplicative noise

?

### n-D Simulations, n=50, d=1:500

#### A) no noise,

$$
rdcorr \approx mdcorr \> HHG
$$

#### B) additive multivariate noise

$$
rdcorr \geq mdcorr \> dcorr,HHG
$$

#### C) multiplicative multivariate noise

??

### non-Euclidean settings

#### A) graphs and covariates

rdcorr works?

#### B) shapes and disease

rdcorr works?

Real Data Results
-----------------

#### A) in real data setting of CxP we have

$$
pval(rdcorr(k\*) < pval(rdcorr) \approx pval(HHG) \approx 0.3 <
pval(mdcorr/dcorr)
$$

#### B) in real data settings of shapes and diseases

?

Discussion
==========

Summary
-------

Related Work
------------

[@Heller12]

 

Extensions
----------

1.  faster (approximate) algorithms

2.  choose optimal k-neighborhood, improve speed and performance

3.  address conditional dependencies

Appendix
========

Proof 1
-------

Proof 2
-------

etc.

Other Misc Notes
================

The draft is named as RDC.pdf in the draft folder.

TibsSimu2Dim is the Matlab code for the first numerical section in RDC.
TibsSimu2Noise is the Matlab code for the second numerical section in RDC.
canoncorr suppress the full-rank warning for the above two codes; but it is not
really used to compare rankdCorr and dCorr.

The pdf and the figures include 1) the 10 tests with respect to dimension from
d=1:1000 without noise (a noisy version looks similar), 2) the 10 test with
respect to Gaussian noise at d=500,

for comparison, the figures further include 3) the 10 test with respect to
Gaussian noise at d=1 (which is similar to the Tibs comments) 4) the 10 test
with respect to Gaussian noise at d=300, which are not included in the pdf
write-up yet.
