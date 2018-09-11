for(i in (1:465)){
  queue = "day"
  mem = "32G"

    com = sprintf("sbatch --mem=%s --partition=%s --wrap=\"R CMD BATCH '--args %s' step5_ase_cnt.R step5_ase_cnt_%s.Rout\"",mem,queue,i,i)
    message(com)
    system(com)
}
