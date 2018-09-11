
library("GenomicFeatures")

root.dir="/lustre/scr/z/h/zhabotyn/R01"
pipe.dir = sprintf("%s/pipe", pipe.dir)
gtf.dir = sprintf("%s/bowtie2_indexes", root.dir)
setwd(pipe.dir)

gtf.file = sprintf("%s/gencode.v21.annotation.gtf")

txdb2 = makeTranscriptDbFromGFF(file=gtf.file, format="gtf",
  exonRankAttributeName="exon_number",
  dataSource=paste("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/",
  "gencode.v21.annotation.gtf.gz",sep=""),
  species="Homo sapiens")

saveDb(txdb2,file=sprintf("%s/Homo_sapiens_GRCh38.sqlite", gtf.dir)

q(save="no")
