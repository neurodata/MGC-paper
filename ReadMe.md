# Local Graph Correlation

some related links:

 - http://emotion.technion.ac.il/~gorfinm/files/science6.pdf - shows that HHG > MIC
 - http://bib.oxfordjournals.org/content/early/2013/08/20/bib.bbt051.short shows again that HHG is the best
 - http://cran.r-project.org/web/packages/HHG/index.html - R package for computing HHG
 - http://arxiv.org/pdf/1308.1559v3.pdf - non-data adaptive partition based test.  seems related to the original idea of partitioning the data up.  we might want to come back to this after publishing the first paper, when we make things zippy by using something like cover trees to partition the data.  actually, they also have a data-derived variant, but it uses a grid, and is gonna be hella ineffective for d>2.
 - http://arxiv.org/abs/1505.02213 - says MICe is "equitable" 
 - http://arxiv.org/abs/1505.02214 - says MICe = HHG > dCov, only consider p=q=1.
 - http://arxiv.org/abs/1410.1503 - fast dCov
 - http://www.exploredata.net/ftp/empirical_ supplement.zip - link to all code used in the Reshef empirical comparison
 - http://projecteuclid.org/euclid.aos/1413810731 - partial dCov
 - http://www.tandfonline.com/doi/abs/10.1080/01621459.2012.695654 - feature screening via dCorr learning
 - https://projecteuclid.org/euclid.aop/1378991840 - dCov in metric spaces, might have useful theory for proving stuff with ranks on metric spaces
 - http://www.sciencedirect.com/science/article/pii/S0047259X13000262 - t-test aot permutation for dCorr
 - http://projecteuclid.org/euclid.aoas/1267453935 - extension to high-D
- On Quantifying Dependence: A Framework for Developing Interpretable Measures
