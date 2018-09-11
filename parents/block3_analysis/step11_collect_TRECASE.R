root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/omnidata", root.dir)
pipe.dir = sprintf("%s/omnipipe", root.dir)
trecase.dir = sprintf("%s/TRECASE_MLE", root.dir)
tre.dir = sprintf("%s/parplus_conv",data.dir)
spe.dir = sprintf("%s/out2e5", pipe.dir)
col.dir = sprintf("%s/par_coll", pipe.dir)
setwd(trecase.dir)

if(!file.exists(col.dir))dir.create(col.dir)
setwd(col.dir)

i = 1
trecase = list.files(path=spe.dir,"TRECASE")
trec = list.files(path=spe.dir,"TREC")
trec = trec[-grep("ASE",trec)]
trecm = matrix(NA,nrow=5e4,ncol=14)
trecasem = matrix(NA,nrow=5e4,ncol=17)

#2e5 children
starti1 = starti2 = 1
for(i in 1:length(trecase)){
  trecasei = read.table(sprintf("%s/%s",spe.dir,trecase[i]),as.is=T,header=T)
  treci = read.table(sprintf("%s/%s",spe.dir,trec[i]),as.is=T,header=T)
  if(i==1){
    colnames(trecm)=colnames(treci)
    colnames(trecasem)=colnames(trecasei)  
  }
  if(nrow(trecasei)>0){
  endi1 = starti1+nrow(trecasei)-1
  trecasem[starti1:endi1,]=as.matrix(trecasei)
  starti1 = endi1+1
  }
  if(nrow(treci)>0){
  endi2 = starti2+nrow(treci)-1
  trecm[starti2:endi2,]=as.matrix(treci)
  starti2 = endi2+1
  }
}
trecasem = trecasem[1:endi1,]
trecm = trecm[1:endi2,]
dim(trecasem)
dim(trecm)

trecase2ind = trecasem
trec2ind = trecm
write.csv(trecasem,"sub_trecase_2e5.csv",row.names=F,quote=F)
write.csv(trecm,"sub_trec_2e5.csv",row.names=F,quote=F)

