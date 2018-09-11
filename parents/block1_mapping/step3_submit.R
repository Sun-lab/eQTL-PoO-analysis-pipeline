for(i in 1:465){
    com = paste("sbatch -n 8  --partition=week --wrap=\'R CMD BATCH \"--args ", i ,"\" step3_qc.R step3_qc_",i,".Rout'",sep="")
    system(com)
}
q("no")

