root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/omnidata", root.dir)
pipe.dir = sprintf("%s/omnipipe", root.dir)
trecase.dir = sprintf("%s/TRECASE_MLE", root.dir)
tre.dir = sprintf("%s/parplus_conv",data.dir)
spe.dir = sprintf("%s/out2e5", pipe.dir)
setwd(trecase.dir)

len = rep(NA,22)
for(i in 1:22){
#i = 1
fli = sprintf("%s/gene_info_chr%s.dat",tre.dir,i,header=T)
fli = read.table(fli)
len[i] = nrow(fli)
}
len


#window 2e5
by = 10
bat = "bat"
for(i in 1:22){
  #i = 1
  ni = ceiling(len[i]/by)
  for(j in 1:ni){
    #j = 31
    starti = (j-1)*by+1
    endi = min(j*by,len[i])
    c(starti,endi)
    com = sprintf("TRECASE_MLE --sfile %s/specification.txt --chr %s:%s-%s",spe.dir,i,starti,endi)
    com2 = sprintf("sbatch --output %s/batch_chr%s_%s.out --partition %s --wrap='%s'",spe.dir,i,j,bat,com)
    system(com2)
  }
}

q("no")
