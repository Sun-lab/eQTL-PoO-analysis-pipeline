#we use shapeit version 2 with hg19 aligned 1000G data
#this processing is done for each chromosome
queue = "bat"                                                                               

root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/omnidata", root.dir)
info.dir = sprintf("%s/info", root.dir)
pipe.dir = sprintf("%s/omnipipe", root.dir)
shape.dir = sprintf("%s/shapeit", root.dir)

geno.dir = sprintf("%s/data_genotype/omnireform", data.dir)
phas.dir = sprintf("%s/phasedR",geno.dir)
refinfo = sprintf("%s/1000GP_Phase3",info.dir)

setwd(shape.dir)

all1000G = "1000GP_Phase3" #common reference prefix
mem = "8G"
for(i in c(1:22)){
  #i = 22;
  if(i %in% 1:8)mem = "16G"
  pref = sprintf("./shapeit -check -T 8 ")
  bedi = sprintf("-B %s/ALL.chr%s.omni_2123_samples_b37_SHAPEIT.20120103.snps.chip_based.haplotypes",geno.dir,i)

  mapr = sprintf("-M %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)
#  inpr = sprintf("--input-ref %s/%s_chr%s_impute.hap.gz",refinfo,all1000G,i)
#  legr = sprintf("%s/%s_chr%s_impute.legend.gz",refinfo,all1000G,i)
  inpr = sprintf("--input-ref %s/%s_chr%s.hap.gz",refinfo,all1000G,i)
  legr = sprintf("%s/%s_chr%s.legend.gz",refinfo,all1000G,i)
  samr = sprintf("%s/%s.sample",refinfo,all1000G)
  #this option is not actually used in the check-step, it will be produced during a check
  #excl = sprintf("--exclude-snp %s/chri/gwas.alignments_%s.snp.strand.exclude",geno.dir,i)
  out = sprintf("--output-log %s/gwas.alignments_%s",geno.dir,i)
  com = sprintf("%s %s %s %s %s %s %s",
                          pref, bedi, mapr, inpr, legr, samr, out)
#  com2 = sprintf("bsub -q %s '%s'",queue,com)
  qout = sprintf("--output=%s/rout/step1_out_alt_%s.out",pipe.dir,i)    
  com2 = sprintf("sbatch %s --partition=%s --mem=%s --wrap='%s'",qout,queue,mem,com)
  message(com)
#  system(com)
  system(com2)
}


q("no")
