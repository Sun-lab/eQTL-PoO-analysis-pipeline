#now when we have exclusion set of snps in step 2 we can run shapeit
queue = "bat"                                                                               
root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/data", root.dir)
info.dir = sprintf("%s/info", root.dir)
pipe.dir = sprintf("%s/pipe", root.dir)
shape.dir = sprintf("%s/shapeit", root.dir)

geno.dir = sprintf("%s/data_genotype/ped90reform", data.dir)
phas.dir = sprintf("%s/phasedR",geno.dir)
if(!file.exists(phas.dir))dir.create(phas.dir)
refinfo = sprintf("%s/1000GP_Phase3",info.dir)

setwd(shape.dir)

all1000G = "1000GP_Phase3" #common reference prefix
mem = "96G"
for(i in c(1:22)){
  #i = 16;
  pref = sprintf("./shapeit")
  pedi = sprintf("--input-ped %s/chri/geno_chr%s",geno.dir,i)
  mapr = sprintf("--input-map %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)
  inpr = sprintf("--input-ref %s/%s_chr%s.hap.gz",refinfo,all1000G,i)
  legr = sprintf("%s/%s_chr%s.legend.gz",refinfo,all1000G,i)
  samr = sprintf("%s/%s.sample",refinfo,all1000G)
  excl = sprintf("--exclude-snp %s/chri/gwas.alignments_%s.snp.strand.exclude",geno.dir,i)
  out = sprintf("-O %s/phasing_chr%s.log",phas.dir,i)
  oth = "--effective-size 20000 --seed 1234567 --thread 8"
  com = sprintf("%s %s %s %s %s %s %s %s %s",
               pref, pedi, mapr, inpr, legr, samr, excl, out, oth)
  qout = sprintf("--output=%s/rout/step4_out_%s.out",pipe.dir,i)    
  com2 = sprintf("sbatch %s --partition=%s --mem=%s --wrap='%s'",qout,queue,mem,com)
  message(com)
  #system(com2)
}


q("no")
