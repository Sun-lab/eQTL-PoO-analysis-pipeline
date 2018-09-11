args = commandArgs(TRUE)
#args = c("1")
indj = as.numeric(args[1])

#we collect the data at the individual level
#also we perform lifting over to to HG38
#ensure that when we write positions to a file R doesn't switch them to scientific
options("scipen"=999,"digits"=4)

root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/omnidata", root.dir)
info.dir = sprintf("%s/info", root.dir)
pipe.dir = sprintf("%s/omnipipe", root.dir)
shape.dir = sprintf("%s/shapeit", root.dir)

geno.dir = sprintf("%s/data_genotype/omnireform", data.dir)
phas.dir = sprintf("%s/phasedR",geno.dir)
refinfo = sprintf("%s/1000GP_Phase3",info.dir)
impu.dir = sprintf("%s/imputed",geno.dir)

out.dir = sprintf("%s/snps_omni",data.dir)
if(!file.exists(out.dir))dir.create(out.dir)
setwd(impu.dir)

lift.dir = sprintf("%s/liftover",root.dir)
setwd(lift.dir)
poss_chr = sprintf("chr%s",c(1:22,"X","Y","M"))

fls = list.files(path=out.dir,pattern="^ind")
summ = matrix(0,nrow=length(fls),ncol=length(poss_chr))
colnames(summ) = poss_chr
sht = matrix(unlist(strsplit(fls,"_")),nrow=2)[2,]
rownames(summ) = sht
unma = summ

j = indj
  chrs = list.files(path=sprintf("%s/%s",out.dir,fls[j]),pattern="chr")
  chrs = chrs[match(sprintf("%s.txt",poss_chr),chrs)]
  chrs = chrs[!is.na(chrs)]
  chrs
  app = F
  i = 1
  for(i in 1:length(chrs)){
    chri = read.table(sprintf("%s/%s/%s",out.dir,fls[j],chrs[i]),as.is=T)
    chri[1:2,]
    if(nrow(chri)!=length(unique(chri[,2])))message("problem with duplicated snp")
    #a quick check whether we've got duplicated snps while combining those imputed values
    chri = aggregate(rep(1,nrow(chri)),by=list(chri[,1],chri[,2],chri[,3],chri[,4],chri[,5], chri[,6]),FUN=sum)[,-7]
    chri = chri[order(chri[,2]),]
    chri0 = aggregate(rep(1,nrow(chri)),by=list(chri[,1],chri[,2]),FUN=sum)
    table(chri0[,3])
    m = (chri[,2] %in% chri0[chri0[,3]>1,2])
    chri = chri[!m,]
    if(nrow(chri)!=length(unique(chri[,2])))message("problem with duplicated snp not resolved")

    #convert to hg38 using liftover 
    #https://mathgen.stats.ox.ac.uk/impute/README_1000G_phase1interim_jun2011.txt
    chri$end = chri[,2]+1
    chri[,2] = chri[,2]
    chri = chri[,c(1,2,7,3,4,5,6)]
    chri[1:2,]
    
    
    inpf = sprintf("%s/tmp_%s.bed", lift.dir, indj)
    inp = sprintf("./liftOver %s/tmp_%s.bed", lift.dir, indj)
    cha = sprintf("%s/chains/hg19ToHg38.over.chain",lift.dir)
    out = sprintf("%s/tmp_hg38_%s.bed",lift.dir, indj)
    unl = sprintf("%s/unm_hg38_%s.bed",lift.dir, indj)
    inp
    cha
    out
    unl
    write.table(chri[,1:3],file=inpf, sep=" ",
                quote=F,row.names=F,col.names=F)
    comm = sprintf("%s %s %s %s",inp, cha,out,unl)
    comm
    system(comm)

    chri2 = read.table(file=out)
    info = file.info(unl)
    if(info$size>0){
      unmi2 = read.table(file=unl)
      to.rm = match(unmi2[,2],chri[,2])
      if(length(to.rm)>0)chri = chri[-to.rm,]
      unma[j,i] = unma[j,i] + nrow(unmi2)
    }
    dim(chri)
    dim(chri2)
    chri[1:5,]
    chri2[1:5,]
    #update positions to hg38
    chri[,2] = chri2[,2]    
    #once liftover is completed we again may have a situation with multiple snps at one position
    #if they don't conflict keep one of them, if they conflict, remove all
    if(nrow(chri)!=length(unique(chri[,2]))){
      to.fix = names(table(chri[,2])[table(chri[,2])>1])
      to.rm = numeric(0)
      for(fixi in to.fix){
        ambig = which(chri[,2]==fixi)
        attempt = aggregate(rep(0,length(ambig)),
        by=list(chri[ambig,1],chri[ambig,2],chri[ambig,5],chri[ambig,6]),FUN=sum)[,-5]
        if(nrow(attempt)==1){
          to.rm = c(to.rm,ambig[-1])
        }else{
          to.rm = c(to.rm,ambig)
        }
      }
      chri=chri[-to.rm,]
      message("removed ",length(to.rm)," ambiguous snps")
    }
    #append to the combined list of snps for this individuals
    write.table(chri[,c(1,2,4,5)],file=sprintf("%s/%s.bed",out.dir,sht[j]),
                quote=F,row.names=F,col.names=F,append=app)
    write.table(chri[,c(1,2,4,5,6,7)],file=sprintf("%s/%s_withid.bed",out.dir,sht[j]),
                quote=F,row.names=F,col.names=F,append=app)
    summ[j,i] = summ[j,i] + nrow(chri)
    app = T    
 }
 message("processed ",j,"'th individual: ",sht[j])
#cleaning up temporary files
file.remove(out)
file.remove(unl)
file.remove(inpf)

q("no")
