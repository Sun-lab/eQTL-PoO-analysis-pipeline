args = commandArgs(trailingOnly=TRUE)
argchri = args[1]
#argchri = 22

root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/omnidata", root.dir)

#par+ data
geno.dirp = sprintf("%s/data_genotype/omnireform", data.dir)
vcf.dirp = sprintf("%s/vcf",geno.dirp)#need

#individual data
pipe.dir = sprintf("%s/pipe", root.dir)
geno.dir = sprintf("%s/data_genotype/ped90reform", data.dir)
phas.dir = sprintf("%s/phasedR",geno.dir)#need
impu.dir = sprintf("%s/imputed",geno.dir)#need
exp.dir = sprintf("%s/snpinfoexp",impu.dir)#need
col.dir = sprintf("%s/collected", pipe.dir)

if(!file.exists(exp.dir))dir.create(exp.dir)

setwd(col.dir)

get_block = function(string,split=" ", block=2){
  unlist(strsplit(string,split=split))[block]
}

chri = 22
#the list of samples for which we will get snps
samples = read.table(sprintf("%s/phasing_chr%s.log.sample",phas.dir,chri),header=T,as.is=T)
father = samples$father[-1]
mother = samples$mother[-1]
samples = samples[-1,1]
#define snps (I'll exclude indels with matching to the base
bases=c("A","C","G","T")


library(qvalue)

trecasem = read.csv("sub_trecase_2e5.csv",as.is=T)
trecm = read.csv("sub_trec_2e5.csv",as.is=T)

trecasem[,"p_perm"] = as.numeric(trecasem[,"p_perm"])
trecasem[,"p_cs.trs"] = as.numeric(trecasem[,"p_cs.trs"])
trecasem[,"minpSNP_Rsq"] = as.numeric(trecasem[,"minpSNP_Rsq"])

trecm[,"p_perm"] = as.numeric(trecm[,"p_perm"])
trecm[,"p_cs.trs"] = as.numeric(trecm[,"p_cs.trs"])
trecm[,"minpSNP_Rsq"] = as.numeric(trecm[,"minpSNP_Rsq"])


dim(trecasem)
dim(trecm)
alpha2 = 0.01

sgngns = union(trecm[,1],trecasem[,1])
sgntraTO = trecm[which(trecm$p_cs.trs<alpha2),1]
sgntraTA = trecasem[which(trecasem$p_cs.trs<alpha2),1]
length(sgntraTA)
length(sgntraTO)

combine = data.frame(genes = sgngns, trn = (sgngns %in% sgntraTA),trnTO = (sgngns %in% sgntraTO),
snp = trecasem$minpSNP_id[match(sgngns,trecasem[,1])],snpPos = trecasem$minpSNP_pos[match(sgngns,trecasem[,1])],
snpTO = trecm$minpSNP_id[match(sgngns,trecm[,1])],snpPos = trecm$minpSNP_pos[match(sgngns,trecm[,1])])

kp1 = (!combine$trn&!combine$trnTO)&!is.na(combine$snp)
kp2 = (!combine$trn&!combine$trnTO)&is.na(combine$snp)
table(kp1)
table(kp2)
table(combine$trn[kp1]+combine$trnTO[kp1])
table(combine$trn[kp2]+combine$trnTO[kp2])
kp3 = (!combine$trn&combine$trnTO)&!is.na(combine$snp)
kp4 = (combine$trn&!combine$trnTO)
table(kp3)
table(kp4)
table(combine$trn[kp3]+combine$trnTO[kp3])
table(combine$trn[kp4]+combine$trnTO[kp4])

write.csv(combine[kp1,],"no_trans_both.csv",quote=F,row.names=F)
write.csv(combine[kp3,],"no_trans_TRECASE.csv",quote=F,row.names=F)
write.csv(combine[kp2|kp4,],"no_trans_TREC.csv",quote=F,row.names=F)

combine$snp[kp1][1:5]
chri = as.numeric(argchri)
nms = read.table(sprintf("%s/SNP_VCF_chr%s.vcf",vcf.dirp, chri),nrow=4,comment.char="@",sep="\t",as.is=T)
nms = nms[4,1]
nms = unlist(strsplit(gsub("#","",nms),split=" "))
nms
if(!file.exists("no_trans_both_ord.csv")){
for(chri in 1:22){
#chri = 22
  snpi = read.table(sprintf("%s/SNP_VCF_chr%s.vcf",vcf.dirp,chri))
  kpf = na.omit(match(combine$snp[kp1],snpi[,3]))
  kpp1 = na.omit(match(combine$snp[kp3],snpi[,3]))
  kpp2 = na.omit(match(combine$snpTO[kp2|kp4],snpi[,3]))  
  message(length(kpf)," ",length(kpp1)," ",length(kpp2))
  
  if(chri == 1){
    snpf = snpi[kpf,1:5]
    snpp1 = snpi[kpp1,1:5]
    snpp2 = snpi[kpp2,1:5]
  }else{
    snpf = rbind(snpf,snpi[kpf,1:5])
    snpp1 = rbind(snpp1,snpi[kpp1,1:5])  
    snpp2 = rbind(snpp2,snpi[kpp2,1:5])  
  }  
  message(chri," ",nrow(snpf)," ",nrow(snpp1)," ",nrow(snpp2))
}
message(sum(kp1)," ",sum(kp3)," ",sum(kp2|kp4))

m = match(snpf[,3],combine$snp[kp1])
table(combine$snp[kp1][m]==as.character(snpf[,3]))
table(combine$snpPos[kp1][m]==as.character(snpf[,2]))


of = order(snpf[,1],snpf[,2])
op1 = order(snpp1[,1],snpp1[,2])
op2 = order(snpp2[,1],snpp2[,2])

snpf = snpf[of,]
snpp1 = snpp1[op1,]
snpp2 = snpp2[op2,]

mf = match(snpf[,3],combine$snp[kp1])
mp1 = match(snpp1[,3],combine$snp[kp3])
mp2 = match(snpp2[,3],combine$snpTO[kp2|kp4])
mf[1:4]
mp1[1:4]
mp2[1:4]

table(as.character(combine$snp[kp1][mf])==snpf[,3])
table(as.character(combine$snp[kp3][mp1])==snpp1[,3])
table(as.character(combine$snpTO[kp2|kp4][mp2])==snpp2[,3])
no_both = combine[kp1,][mf,]
no_TRECASE = combine[kp3,][mp1,]
no_TREC = combine[kp2|kp4,][mp2,]

snpf$Gene_name = no_both$genes
snpp1$Gene_name = no_TRECASE$genes
snpp2$Gene_name = no_TREC$genes
write.csv(cbind(no_both,snpf),"no_trans_both_ord.csv",quote=F,row.names=F)
write.csv(cbind(no_TRECASE,snpp1),"no_trans_TRECASE_ord.csv",quote=F,row.names=F)
write.csv(cbind(no_TREC,snpp2),"no_trans_TREC_ord.csv",quote=F,row.names=F)

conf = no_both
con1 = no_TRECASE
con2 = no_TREC
}else{
  conf = read.csv("no_trans_both_ord.csv",as.is=T)
  con1 = read.csv("no_trans_TRECASE_ord.csv",as.is=T)
  con2 = read.csv("no_trans_TREC_ord.csv",as.is=T)

  snpf =  conf[,-(1:7)]
  snpp1 = con1[,-(1:7)]
  snpp2 = con2[,-(1:7)]
}

snpf[,3] = as.character(snpf[,3])
snpp1[,3] = as.character(snpp1[,3])
snpp2[,3] = as.character(snpp2[,3])

chri = as.numeric(argchri)
snpip = c(snpf[snpf[,1]==sprintf("chr%s",chri),2],snpp1[snpp1[,1]==sprintf("chr%s",chri),2],snpp2[snpp2[,1]==sprintf("chr%s",chri),2])
snpi = rbind(snpf[snpf[,1]==sprintf("chr%s",chri),],snpp1[snpp1[,1]==sprintf("chr%s",chri),],snpp2[snpp2[,1]==sprintf("chr%s",chri),])
snpi[,3] = gsub("chr","",snpi[,3])
snpish2 = sapply(snpi[,3],get_block,split=":",block=2)
snpish = sapply(snpi[,3],get_block,split=":",block=1)
snpish[snpish==chri] = snpish2[snpish==chri]
snpish2[is.na(snpish2)] = snpish[is.na(snpish2)]
names(snpish2)=names(snpish)=NULL
snpish

ngen = nrow(snpi)
geni = rep(NA,ngen)

#get all the appropriate chunks for the given chromosome of trios
#here we use direct results of ped/map format
fls = list.files(path=impu.dir,pattern=sprintf("chr%s_",chri))
to.rm = union(union(union(grep(".txt",fls),grep("info",fls)),grep("warnings",fls)),grep("summary",fls))
fls = fls[-to.rm]
hps = fls[grep("_haps",fls)]
alp = fls[grep("_allele_probs",fls)]
fls = setdiff(fls,union(hps,alp))
fls
hps
alp

i = 1
nfls = length(alp)
for(i in 1:nfls){
  alpi = readLines(sprintf("%s/%s",impu.dir,hps[i]))
  alpiid = sapply(alpi, get_block)
  alpips = sapply(alpi, get_block, block=3)
  names(alpiid) = NULL
  alpiidsh2 = sapply(alpiid,get_block,split=":",block=2)
  alpiidsh = sapply(alpiid,get_block,split=":",block=1)
  alpiidsh[alpiidsh==chri] = alpiidsh2[alpiidsh==chri]
  alpiidsh2[is.na(alpiidsh2)] = alpiidsh[is.na(alpiidsh2)]
  names(alpiidsh2)=names(alpiidsh)=NULL
  alpiidsh[1:5]
  alpiidsh2[1:5]
  
  m1 = match(snpish,alpiidsh)
  m2 = match(snpish,alpiidsh2)
  m3 = match(snpish2,alpiidsh)
  m4 = match(snpish2,alpiidsh2)
  m5 = match(snpish2,alpips)

  geni[!is.na(m5)] = alpi[na.omit(m5)]  
  geni[!is.na(m4)] = alpi[na.omit(m4)]  
  geni[!is.na(m3)] = alpi[na.omit(m3)]
  geni[!is.na(m2)] = alpi[na.omit(m2)]
  geni[!is.na(m1)] = alpi[na.omit(m1)]
  message(i, " out of ", nfls, " : ", sum(is.na(geni)), " out of ", ngen)
}
table(is.na(geni))
#snpdir = "snpinfoexp"
write.table(cbind(snpi,geni),sprintf("%s/expand_chr%s.txt",exp.dir,chri),row.names=F,col.names=F,quote=F)

q("no")