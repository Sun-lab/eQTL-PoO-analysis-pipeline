#we use shapeit version 2 with hg19 aligned 1000G data
#this processing is done for each chromosome
#at this stage we check if the snps labeled as "strand issue" 
#can be fixed by flipping them with plink
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

checki = matrix(NA,nrow=22,ncol=2)
rownames(checki)=sprintf("chr%s",1:22)
#as with each shapeit we do it on chromosome level
all1000G = "ALL_1000G_phase1integrated_v3" #common reference prefix
for(i in c(1:22)){
  #i = 22;

  #flip the snps of interest
  con = file(sprintf("%s/chri/gwas.alignments_%s.snp.strand",geno.dir,i))
  lns = readLines(con)
  lns = lns[substr(lns,1,6)=="Strand"]
  strand = matrix(unlist(strsplit(lns,split="\t")),nrow=11)[4,]
  close(con)
  write.table(unique(strand),"tmp_strand.txt",row.names=F,col.names=F,quote=F)
#  bedi = sprintf("-B %s/ALL.chr%s.omni_2123_samples_b37_SHAPEIT.20120103.snps.chip_based.haplotypes",geno.dir,i)

#  mapr = sprintf("-M %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)

  file.copy(from = sprintf(s/ALL.chr%s.omni_2123_samples_b37_SHAPEIT.20120103.snps.chip_based.haplotypes.bed",geno.dir,i), to = "tmp.bed")
  file.copy(from = sprintfs/genetic_map_chr%s_combined_b37.txt",geno.dir,i), to = "tmp.map")
  system("plink --file tmp --flip tmp_strand.txt  --recode")
  cat("plink done ")
  
  con = file(sprintf("tmp.bed",geno.dir,i))
  tmp = readLines(con)
  close(con)
  cat("\t read preconv ")

  con = file(sprintf("plink.bed",geno.dir,i))
  tmp0 = readLines(con)
  close(con)
  cat("\t read postconv: ")
  #check that the plink really flipped those snps
  #we see that plink does flip those strands
  message(!all(unlist(tmp)==unlist(tmp0)))
  
  
  #recheck using shapeit on flipped data
  pref = sprintf("./shapeit -check")
  bedi = sprintf("-B plink")
  mapr = sprintf("-M %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)
  inpr = sprintf("--input-ref %s/%s_chr%s_impute.hap.gz",refinfo,all1000G,i)
  legr = sprintf("%s/%s_chr%s_impute.legend.gz",refinfo,all1000G,i)
  samr = sprintf("%s/%s.sample",refinfo,all1000G)
#  excl = sprintf("--exclude-snp %s/chri/gwas.alignments_%s.snp.strand.exclude",geno.dir,i)
  out = sprintf("--output-log tmp")
  com = sprintf("%s %s %s %s %s %s %s",
                          pref, bedi, mapr, inpr, legr, samr, out)
  message(com) 
  system(com)
  
  con = file(sprintf("tmp.snp.strand",geno.dir,i))
  lns2 = readLines(con)
  lns2 = lns[substr(lns,1,6)=="Strand"]
  strand2 = matrix(unlist(strsplit(lns2,split="\t")),nrow=11)[4,]
  close(con)
  
  #count how many snps are in the one exclusion set, but not in the other
  checki[i,1] = length(setdiff(unique(strand),unique(strand2)))
  checki[i,2] = length(setdiff(unique(strand2),unique(strand)))
  message("processed ",i,"'th file")
}

#remove temporary files
system("rm -f tmp*")
system("rm -f plink*")

#we don't see any snps resolved this way, file contains zeroes
write.table(checki,sprintf("%s/step2_check.txt",pipe.dir),
            row.names=T,col.names=F,quote=F)

q("no")
