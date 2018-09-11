root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/data",root.dir)
geno.dir = sprintf("%s/data_genotype/ped90reform", data.dir)
vcf.dir = sprintf("%s/vcf",geno.dir)
setwd(vcf.dir)
chr.dir = sprintf("%s/chrs",vcf.dir)
if(!file.exists(chr.dir))dir.create(chr.dir)
vcf.dir
vcfs = list.files()
vcfs = vcfs[-grep("chrs",vcfs)]
templ = vcfs[grep("templ",vcfs)]
vcfs = vcfs[-grep("templ",vcfs)]

get_block = function(str,split="_",block=1){
  unlist(strsplit(str,split=split))[block]
}
for(i in 1:22){
  chrid = sapply(vcfs,get_block,block=3)
  vcfi = vcfs[chrid==sprintf("chr%s",i)]
  o = order(as.numeric(sapply(vcfi,get_block,block=4)))
  
  inp = paste(vcfi[o],collapse = " ")
  out = sprintf("> %s/SNP_VCF_chr%s.vcf",chr.dir,i)
  com = paste("cat",templ,inp,out,sep=" ")
  com
  system(com)
  message(i)
}
