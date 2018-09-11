for(i in 1:465){
    com = paste("sbatch -n 8  --partition=week --wrap=\'R CMD BATCH \"--args ", i ,"\" step4_check_and_filter_liberal.R step4_check_and_filter_liberal_",i,".Rout'",sep="")
    system(com)
}
q("no")

