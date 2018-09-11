if(!file.exists("rout"))dir.create("rout")
for(i in (1:22)){
#i = 1
  mem = "16G"
  queue = "p01_bat"
  com = sprintf("R312 CMD BATCH '--args %s' step6_get_snps_upd.R rout/step6_get_snps_upd_%s.Rout", i, i)
  qout = sprintf("--output=rout/outsh_%s.out", i)
  com2 = sprintf("sbatch %s --partition=%s --mem=%s --wrap=\"%s\"",qout,queue,mem,com)
  message(com)
  system(com2)     
}

q("no")
