#getting a chromosome level summary of how many snps are to be removed 
#due to missing values or to strand issue

root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/data", root.dir)
info.dir = sprintf("%s/info", root.dir)
pipe.dir = sprintf("%s/pipe", root.dir)

geno.dir = sprintf("%s/data_genotype/ped90reform", data.dir)
phas.dir = sprintf("%s/phasedR",geno.dir)
refinfo = sprintf("%s/1000GP_Phase3",info.dir)
setwd(pipe.dir)

  reform.dir = sprintf("%s/chri%s",geno.dir)
  
  chrs = 1:22
  summ = matrix(NA,nrow=length(chrs),ncol=4)
  colnames(summ) = c("total","not.comp","strand","missing")
  rownames(summ) = paste("chr",chrs,sep="")
  for(i in chrs){
    #i = 16
    geno = read.table(sprintf("%s/geno_chr%s.map",reform.dir,i),as.is=T)
    summ[i,1] = nrow(geno)
    
    dat = read.table(sprintf("%s/gwas.alignments_%s.snp.strand",
                     reform.dir,i),sep=",",header=F,as.is=T)
    dat = dat[-1,]
    summ[i,2] = length(dat)
    tbl = table(matrix(unlist(strsplit(dat,"\t")),ncol=length(dat))[1,])
    summ[i,3] = tbl[["Strand"]]
    summ[i,4] = tbl[["Missing"]]
    message(i)
  }
  write.csv(summ, "step3_snp_summary.csv", quote=F)
  message(suff)

q("no")
