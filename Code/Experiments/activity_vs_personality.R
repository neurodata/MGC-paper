# youngser park (2013), called this function "CGP" for "connectome x genome x personality"
# this code pre-processed the data for the activity vs personality test

library(XML)
library(R.matlab)

dname
<-
"http://www.cis.jhu.edu/~parky/CGP/Data/time_series_cc200_20130402/"
urlfiles
<-
readHTMLTable(dname,skip.rows=1:2)[[1]]$V2
urlfiles
<-
urlfiles[!is.na(urlfiles)]
files
<-
paste0(dname,urlfiles)
(n
<-
length(files))
# 42
require(flexmix)
coh
<-
ts
<-
spect
<-
spect_norm
<-
NULL
for
(i
in
1:length(files))
{
ts[[i]]
<-
readMat(files[i])[[1]]
# 194 x 197
spect[[i]]
<-
spectrum(t(ts[[i]]))$spec[6:40,]
# 100 (freq) x 194
 => (bandpass) 35 x 194
spect_norm[[i]]
<-
scale(spect[[i]],center=FALSE,scale=colSums(spect[[i]]))
# normalization
coh[[i]]
<-
KLdiv(spect_norm[[i]])
# Delta: 194 x 194
}
