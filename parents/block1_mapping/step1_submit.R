for(i in 1:465){
    com = paste("sbatch -n 8  --partition=week --wrap=\'R CMD BATCH \"--args ", i ,"\" step1_trim_filter.R step1_trim_filter_",i,".Rout'",sep="")
    system(com)
}
q("no")

