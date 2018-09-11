#convert samples
setwd("/home/groups/projects/Sun_RNA_seq/data_genotype/omnireform/")
for(chri in 1:22){
   bnm = sprintf("ALL.chr%s.omni_2123_samples_b37_SHAPEIT.20120103.snps.chip_based.haplotypes", chri)

  com = sprintf("plink --vcf gz/%s.vcf.gz --make-bed --out %s",bnm, bnm)
  system(com)
}


q("no")
