---
title: Local Independence Testing
author: sigma(Cencheng Shen, cep, mauro, jovo)
---

# Abstract

Understanding and discovering dependence between multiple properties or measurements of our world is a fundamental task not just in science, but also policy, commerce, and other domains. In the past hundred years, people have developed many different measures of dependence that can be applied in a wide variety of settings.  An ideal dependence measure would have the following properties. (1) Strong theoretical support, guaranteeing rejecting independence no matter what the dependence structure is. (2) Strong empirical support on a wide variety of low- and high-dimensional simulation settings. (3) Provides insight into the local scale in which dependency is strongest. (4) Detects dependence when it exists, and fails to detect dependence when it does not exist, on real data. No existing test satisfies all of these properties. We develop a novel dependence statistic and test called "Local Graph Dependence" that does. We can therefore use this test in a variety of settings in which previous tests failed to detect signal or provide insight.

# Outline

GAP: It is easy to collect many different measurements/views/features on different datasets, such as height & weight, genome & personality, structural & functional connectome. 
Lots of work on independence tests, including tests that work on linear data, nonlinear data, and high-dimensional data.
However, there is no test that has all the following properties:
A. strong theoretical support
B. empirically outperforms others on low-dimensional settings
C. empirically outperforms others on high-dimensional settings
D. provides insight into local structure of data.
E. empirically outperforms other stuff on real data


CLAIM: we have a test  that satisfies A-E.

Figure 1: On 20 different simulation settings, using the exact parameters proposed in previous papers, our method empirically achieves as high or higher power than competing approaches for nearly all sample sizes on nearly all problems (CLAIM B).

Figure 2: Quantitatively evaluating the performance of the various algorithms on the 20 benchmarks into a single number, it is clear that our method is overall much better, both (A) for a given sample size for 1D,  (B) across all sample sizes for 1D (CLAIM B), (C) for all dimensions for fixed sample size, and all dimensions and all sample sizes (CLAIM C).

Figure 3: On the "same" 20 different simulation settings, but now increasing the dimensionality from 1 up to 1000, our method empirically achieves as high or higher power than competing approaches for nearly all dimensions on nearly all problems (CLAIM C).

Figure 4: For each of the 20 different simulations, our method (and only our method) provides an estimate of the local scales that encodes the dependence structure.

Figure 5:  On 4 different real data examples, our methods achieves the expected outcome, whereas the others do not.  This means obtaining a significant p-value on 3 settings where we expect dependence should exist (but had never previously been detected), and a 4th where we expect it should not (but previous tests indeed detected a signal). (CLAIM D)

Supp Figure 1: workflow
Supp Figure 2: Example of 100 samples from the 20 different functions (1D variants) for visualization purposes.

Each figure has a single declarative sentence with a verb that establishes what that figure offers to the world.

# Results Subsection Titles:

## Theoretical Claims

## 1D Simulation Results

## High-dimensional Simulation Results

## Local Structure

## Benchmark Data



# Introduction

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

# Results


## Theoretical Results


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

[@Heller2012, @Josee2014]

 

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

Folder '1' is for CorrIndTest, including 20 simulation types for increasing sample size at fixed dimension m. It compares rank dcorr with knn, dcorr, mdcorr, and HHG. It also plots the performance profile at n=40 for all simulations (m is different though), clearly rdcorr is the best.

Folder '1 with neighbor' is based on the same results, but plots the power with respect to increasing neighborhood at fixed n and fixed dimension. We observe that for linear datasets the optimal k=n-1, and for nonlinear datasets the optimal k<n-1.

Folder '2' and '2 with neighbor repeats the same plots', but add knn to dcorr and mdcorr as well. We can see mdcorr is comparable to rdcorr this time, both are better than dcorr and HHG. Two performance profiles are n=40 and n=n_{max}/2 are provided. In this sense, if we cannot work out the consistency of rdcorr, mdcorr with knn is still a method with theoretical consistency and huge performance improvement; although the emphasize of the paper may be different if we chose to go that route.

For RealData1, it includes the p-value for Brain CxP, and Shape vs Disease. Clearly rdcorr is the superior method.

For RealData2, it considers dcorr and mdcorr with knn again. All three methods perform quite well now, although rdcorr is still the only method of significance for brain CxP. 
