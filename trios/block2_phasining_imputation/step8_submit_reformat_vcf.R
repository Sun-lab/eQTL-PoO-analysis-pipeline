if(!file.exists("rout"))dir.create(rout)
for(i in (1:90)){
  queue = "week"
  mem = "8G"
  com = sprintf("R CMD BATCH '--args %s' step8_reformat_vcf.R rout/step8_reformat_vcf_%s.Rout",i,i)
  qout = sprintf("--output=rout/outsh8_%s.out", i)
  com2 = sprintf("sbatch %s --partition=%s --mem=%s --wrap=\"%s\"",qout,queue,mem,com)

  message(com)
  system(com)
}

q("no")
