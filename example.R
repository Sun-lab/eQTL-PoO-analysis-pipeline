percase = 0.1
dblcnt = 0.2
mn = 100
b0 = 0;b1 = 0;th = .5;dv=4;niter = 100;betas = c(3,.2,.05,.5);ss=2

#rho = 0
set.seed(12345)
workdir = "/lustre/scr/z/h/zhabotyn/R01/2018_08_15"
setwd(workdir)
library(rcppreqtl, lib.loc="/nas02/home/z/h/zhabotyn/research/package9/")
library(VGAM)
library(MASS)
#source("/lustre/scr/z/h/zhabotyn/R01/2018_08_15/simu.R")
#source("/lustre/scr/z/h/zhabotyn/R01/2018_08_15/fit.R")
ls("package:rcppreqtl")


#
#setup
#
phiNB = th
phiBB = th/dv
rho = phiBB/(1+phiBB)
c(phiNB, phiBB, rho, th)
dep = makeXmatr(ss)

dat = simu4(num=niter, Xmatr=dep$Xmatr, haplotype=dep$thp, totmean=mn, percase=percase, dblcnt=dblcnt, phiNB=phiNB, phiBB=phiBB, b0=b0, b1=b1, betas=betas)

fullest = fit(subset=NULL, data=dat, traceit=FALSE)
commest = fitsh(subset=NULL, data=dat, traceit=FALSE)

dat0 = dat
dat0$asn=dat0$asnA
dat0$asnp=dat0$asnpA
dat0$haplotype=dat0$haplotypeA
fullest4 = fit(subset=NULL, data=dat0, traceit=FALSE)
commest4 = fitsh(subset=NULL, data=dat0, traceit=FALSE)

mean(fullest$full[,"pval_b0"]<.05)
mean(commest$full[,"pval_b0"]<.05)
mean(fullest4$full[,"pval_b0"]<.05)
mean(commest4$full[,"pval_b0"]<.05)




