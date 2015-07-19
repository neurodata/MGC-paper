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

where f includes 20 types of functions, x is $d \times n$, C is $1\times d$
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

2.  dcorr is O(N\^2), HHG is O(N\^2 log N), rdcorr is O(N\^2 log N), 

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

4.  incorporate  "Modified distance covariance statistics" from Szekeley et al, 2013 to work in metric data, which might be better?


Appendix
========

Proof 1
-------

Proof 2
-------

etc.

Cencheng's Temporary Notes
================

The draft is named as RDC.pdf in the draft folder, which is not yet updated.

The data folder includes:

CorrIndTest is the main Matlab code for independence test at fixed dimension and increasing sample size. It returns the testing powers for dCorr, rank dCorr, modified dCorr, and HHG, which are implemented separately.

CorrPermTest is the main Matlab code for permutation test at fixed dimension and increasing sample size. It returns the p-value. (The figures are not included yet)

CorrIndTestDim is the Matlab code for independence test at fixed sample size and increasing dimension.

CorrIndTestNoise is the Matlab code for independence test at fixed sample size & dimension with increasing noise.

The subfolders in the data folder contain the figures for the respective tests. CorrPermTest & CorrIndTestNoise plots will be included later.


Numerical Summary:

The code includes 20 simulation types, including 10 similar to Tibshirani comment 2011 (some of them are in high-dimensional setting), 7 from Wikipedia correlation figure/Heller 2012 Table 3, and 3 from Szeley 2007 example 1-3. 
All numerical plots are very comparable to existing plots from the above papers, though the exact experiment set-up is a bit different.

In our numerical experiments, rank dCorr shows robust performance against high-dimensional data and nonlinear data, which usually attains close to optimal powers among all 4 statistics, if not the best. 

For rank dCorr, Type 13/20 are the only two cases where rank dCorr appears significantly lower than the best statistic; but rank dCorr is still much better than dCorr/modified dCorr in those two cases

HHG appears to work well for nonlinear data, but it is usually the worst statistic for high-dimensional data. Examples include type 1-3, 8-9, as well as type 18, which is a simple joint normal distribution.

dCorr and modified dCorr usually perform well for high-dimensional data and linear relationship, but not so much when nonlinearity comes into play. They perform significantly worse than rank dCorr and HHG in type 4-7, type 10-16, type 19-20.
