queue = "week"
for(i in c(1:22)){
  batchout = sprintf("step12_produce_cis_list_%s.Bout",i)
  rcom = sprintf("R CMD BATCH '--args %s' step12_produce_cis_list.R step12_produce_cis_list_%s.Rout",i,i)
  com = sprintf("sbatch --output=%s --mem 32G --partition %s --wrap=\"%s\"", batchout,queue,rcom)
  message(com)
  system(com)
}

q("no")
