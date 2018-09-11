if(!file.exists("rout"))dir.create("rout")
for(i in (1:90)){
#i = 1
  mem = "16G"
  queue = "bat"
  com = sprintf("R CMD BATCH '--args %s' step7_collect_snps.R rout/step7_collect_snps_%s.Rout", i, i)
  qout = sprintf("--output=rout/outsh_%s.out", i)
  com2 = sprintf("sbatch %s --partition=%s --mem=%s --wrap=\"%s\"",qout,queue,mem,com)
  message(com)
  system(com2)     
}

q("no")
