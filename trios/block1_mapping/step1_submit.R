for(k in 1:5){
    com = paste("sbatch -n 8  --partition=week --wrap=\'R CMD BATCH \"--args ",k,"\" step1_trim_filter.R step1_trim_filter_",k,".Rout'",sep="")
    system(com)
}
q("no")

