root.dir = "/lustre/scr/z/h/zhabotyn/R01"
pipe.dir = sprintf("%s/pipe", root.dir)

data.dir = sprintf("%s/data",root.dir)
#geno.dir = sprintf("%s/data_genotype/ped90reform", data.dir)
geno.dir = sprintf("%s/data_genotype/ped90", root.dir)
phas.dir = sprintf("%s/phasedR",geno.dir)
impu.dir = sprintf("%s/imputed",geno.dir)
exp.dir = sprintf("%s/snpinfoexp",impu.dir)
cnt.dir = sprintf("%s/ind_new/prep",data.dir)


colp.dir = sprintf("%s/collected", pipe.dir)
setwd(colp.dir)


conf = read.csv("no_trans_both_ord.csv",as.is=T)
con1 = read.csv("no_trans_TRECASE_ord.csv",as.is=T)
con2 = read.csv("no_trans_TREC_ord.csv",as.is=T)

setwd(pipe.dir)
snpf =  conf[,-(1:7)]
snpp1 = con1[,-(1:7)]
snpp2 = con2[,-(1:7)]
conf[1:5,]
snpf[1:5,]

cons = rbind(conf,con1,con2)
snp = rbind(snpf,snpp1,snpp2)
snp$type = c(rep(1,nrow(snpf)),rep(2,nrow(snpp1)),rep(3,nrow(snpp2)))
snp[,1] = gsub("chr","",snp[,1])

snp_sid = read.table(sprintf("%s/phasing_chr22.log.sample",phas.dir),header=T,as.is=T)
snp_sid = snp_sid[-(1),]
for(chri in 1:22){
  snpi = read.table(sprintf("%s/expand_chr%s.txt",exp.dir,chri),as.is=T)
  snpi[,2] = chri
  if(chri==1){snps = snpi}
  else{snps = rbind(snps,snpi)}
  message(chri)
}
dim(snps)
dim(snp)
dim(cons)
get_block = function(string,split=" ", block=2){
  unlist(strsplit(string,split=split))[block]
}
snp$posold = sapply(snp[,3],get_block,split=":")
sht = sapply(snp[is.na(snp$posold),3],get_block, split=":", block=1)
m = match(sht,sapply(snps[,1],get_block,split=":", block=1))
snp$posold[is.na(snp$posold)] = snps[m,4]

#o = order(as.numeric(snps[,2]),as.numeric(snps[,9]))
o = order(as.character(snps[,6]))
snps = snps[o,]
#o = order(as.numeric(snp[,1]),as.numeric(snp[,7]))
o = order(as.character(snp[,6]))
snp = snp[o,]
dim(snps)
dim(snp)
table(snps[,9]==snp[,8])
table(snps[,6]==snp[,6])
cbind(snps[,c(2,6,9)],snp[,c(1,7,6)])[1:5,]



#easy to see that first haplotype of a child = first haplotype of father and second haplotype of a child = first haplotype in mother:
i = 2
offset = 11
for(i in 1:30){
  fathcol = offset + 6*(i-1) + 1#12 & 13  then 18 & 19
  mothcol = offset + 6*(i-1) + 3#14 & 15 then 20 & 21
  chifcol = offset + 6*(i-1) + 5#16 then 22
  chimcol = offset + 6*(i-1) + 6#17 then 23
  
  mf = snps[,fathcol] == snps[,chifcol]
  mm = snps[,mothcol] == snps[,chimcol]
  if(i == 1){
    m = !(mf&mm)
  }else{
    m = cbind(m,!(mf&mm))
  }
  message("ind ", i, " mismatched father: ",sum(!mf)," mother: ",sum(!mm)," overall: ",sum(!(mf&mm)), " out of: ",length(mf))
}
table(rowSums(m))
#14 snps that had contradictory haplotype information in 1 out of 30 individuals, 
# snp that had contradictory haplotype information in 2 individuals 
#1 with contradictory in 3 individuals
apply(m,2,which)
#considerations show that these are unlikely due to recombination (since nearby snps before and after satisfy the same paternal/maternal assumption, so we rather may consider that this question is due to snp error (will keep them though since don't know where the error occured - in parents or in child and we care only of a child)

#the original snps were produced to reflect father (hap1) and mother (hap2) for allele-specific counts
#also in the snp list of interest we observe consistent story first father is listed, then mother is listed 
#and than father/mother haplotypes for individual


#prb[thp==0] = b1                # both genotype 0, first paternal
#prb[thp==1] = b0 + b1           # this is formulated in terms of gi being BA (and first being paternal)
#prb[thp==2] = -b0 + b1          # this is formulated in terms of gi being AB (and first being paternal)
#prb[thp==3] = b1                # both genotype 1, first paternal

#general procedure is to look at the SNP for a particular individuum
#if we observe 00 or 11 - don't touch counts and code them as 0 and 3
#if we observe 10 - don't touch counts and code them as 1  (BA and first paternal)
#if we observe 01 - don't touch and code them as 2 (AB and first paternal)

snps$type=snp$type
dat_sid = read.table(sprintf("%s/samples.dat",cnt.dir),as.is=T,header=T)
for(chri in 1:22){
  info = read.table(sprintf("%s/gene_info_chr%s.dat",cnt.dir,chri),header=T,as.is=T)
  as1 = read.table(sprintf("%s/gene_a1_chr%s.dat",cnt.dir,chri),header=F,as.is=T)
  as2 = read.table(sprintf("%s/gene_a2_chr%s.dat",cnt.dir,chri),header=F,as.is=T)
  tot = read.table(sprintf("%s/gene_Tcount_chr%s.dat",cnt.dir,chri),header=F,as.is=T)
  kp = !is.na(match(info[,1],snps[,6]))
  if(chri == 1){
    infoi = info[kp,]
    as1i = as1[kp,]
    as2i = as2[kp,]
    toti = tot[kp,]
  }else{
    infoi = rbind(infoi,info[kp,])
    as1i = rbind(as1i,as1[kp,])
    as2i = rbind(as2i,as2[kp,])
    toti = rbind(toti,tot[kp,])
  }
  message(chri)
}
dim(infoi)
dim(snps)
#out of 13154 parent+ genes we see enough counts in 12386 children
#get the corresponding SNPs
#for now skip to VCF from VCF conversion
m = match(infoi[,1],snps[,6])
pat = snps[m,offset + 6*(0:29) + 5]
mat = snps[m,offset + 6*(0:29) + 6]
#write.table(snp[m,], sprintf("%s/subsnpinfo.txt",exp.dir),row.names=F,col.names=T,quote=F)

indx = matrix(0, nrow=nrow(pat), ncol=ncol(pat))
colnames(indx) = snp_sid[(1:30)*3,1]
rownames(indx) = snps[m,6]
indx[pat==1&mat==0] = 1
indx[pat==0&mat==1] = 2
indx[pat==1&mat==1] = 3

m = match(dat_sid[,1],colnames(indx))
indx = indx[,m]
all(colnames(indx)==dat_sid[,1])
#write.table(indx, sprintf("%s/index.txt",exp.dir),row.names=T,col.names=T,quote=F)



indx = read.table(sprintf("%s/index.txt",exp.dir),as.is=T)
indx = as.matrix(indx)

trcm = as.matrix(toti)
asnm = as.matrix(as1i + as2i)
asnpm = as.matrix(as1i)

get_numcl = function(vec){
  length(unique(vec))
}
nclas = apply(indx,1,get_numcl)
table(nclas)
which(nclas==1)


#fit the model
#detach("package:rcppreqtl")
library(rcppreqtl, lib.loc="~/research/package9")
library(VGAM)
library(MASS)
niter = nrow(trcm)

Xmatr0 = cbind(log10(colSums(trcm)),diag(1,3)%x%matrix(1,10))
dat_sid$int = 1
Xmatr = as.matrix(dat_sid[,c(10,2,8,9)])
#write.table(Xmatr, "Xmatr.csv", quote=F, row.names=F, col.names=F, sep=",")
nbeta = ncol(Xmatr)

rownames(trcm) = infoi[,1]
rownames(asnm) = infoi[,1]    
rownames(asnpm) = infoi[,1]    

#formulate dep and dat from real data - TReC & AS counts
#simulation example
dat = list(haplotype=indx, trc=trcm, asn=asnm, asnp=asnpm, X=Xmatr, haplotypeA=NULL, asnA=NULL, asnpA=NULL, params=NULL, settings=NULL)

start = proc.time()
fullest = fit(subset=1:niter, data=dat, traceit=FALSE)
end = proc.time()
end-start
names(fullest)
fullest$full[1:4,]
(end-start)/niter

write.csv(fullest$full, "full_model.csv")
write.csv(fullest$testadd, "short_testadd.csv")
write.csv(fullest$testpoo, "short_testpoo.csv")

#simulation example
dep = makeXmatr(ss)
mn = 100
percase = .1
dblcnt = .1
phiNB = 1
phiBB = .5
b0 = 1
b1 = 1
betas = c(1, .5, .25, .4)
dat = simu4(num=niter, Xmatr=dep$Xmatr, haplotype=dep$thp, totmean=mn, percase=percase, dblcnt=dblcnt, phiNB=phiNB, phiBB=phiBB, b0=b0, b1=b1, betas=betas)
dat2 = list(haplotype=dat$haplotype, trc=dat$trc, asn=dat$asn, asnp=dat$asnp, haplotypeA=NULL, asnA=NULL, asnpA=NULL, X=dep$Xmatr, params=NULL, settings=NULL)


#simulation example
start = proc.time()
fullest = fit(subset=NULL, data=dat, traceit=FALSE)
end = proc.time()
end-start
names(fullest)
fullest$full[1:4,]

q("no")