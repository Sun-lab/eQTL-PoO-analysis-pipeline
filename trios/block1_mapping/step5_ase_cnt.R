args = commandArgs(TRUE)
args
#args=c("31")
#this processing is done for each bam file individually
i = as.numeric(args[1])

#note, before running this file run in killdevil
#module add samtools
root.dir = "/lustre/scr/z/h/zhabotyn/R01"
bow2.dir = sprintf("%s/bowtie2_indexes", root.dir)
snp.dir = sprintf("%s/anno/snps",root.dir)
setwd(root.dir)
data.dir = sprintf("%s/data",root.dir)
count.dir = sprintf("%s/cntu38",data.dir)
if(!file.exists(count.dir))dir.create(count.dir)

batches = c("11_03_14","12_19_14","12_29_14","01_12_15","01_22_15")
batches = sprintf("%s/%s",data.dir,batches)

#create a list of paths and files required for each i'th file
outpath = inpath = fls = character(0)
for(j in 1:length(batches)){
  bam.dir = sprintf("%s/tophat_output_liberal/bam_liberal",batches[j])
  filt.dir = sprintf("%s/tophat_output_liberal/ufilt",batches[j])
  if(!file.exists(filt.dir))dir.create(filt.dir)

  flj = list.files(path=bam.dir,pattern="bam")
  flj = flj[flj!="tmp.bam"]
  to.rm = grep("bai",flj)
  to.rm =union(grep("tmp",flj),to.rm)
  if(length(to.rm)>0)flj = flj[-to.rm]
  message(length(flj))

  inpath = c(inpath,rep(bam.dir,length(flj)))
  outpath = c(outpath,rep(filt.dir,length(flj)))
  fls = c(fls,flj)
}
inpath[i]
outpath[i]
fls[i]
outfls = unlist(strsplit(fls,"_filtered.bam"))
 
#get shorter id in order to match snps produced in step 7
#here I assume snps will have the id in the name, otherwise need to supply a connection
get_block = function(str,block,split="_"){
  unlist(strsplit(str,split=split))[block]
}
shtid = sapply(fls,get_block,split="-",block=2)
shtid = sapply(shtid,get_block,split="_",block=1)
shtid = sapply(shtid,get_block,split="Gm",block=2)
table(shtid)

snplists = list.files(snp.dir)
snpid = sapply(snplists,get_block,split=".bed",block=1)
snpid = sapply(snpid,get_block,split="NA",block=2)

matchsnp = match(shtid,snpid)

#perform filtering and splitting into allele specific counts
library(asSeq)
library(isoform)

fli = sprintf("%s/%s",inpath[i],fls[i])
outi = sprintf("%s/%s",outpath[i],outfls[i])
louti = sprintf("%s_sorted_by_name_uniq_filtered.bam",outi)
outi
if(!file.exists(louti))prepareBAM(fli, outi,min.avgQ = 30, min.mapQ = 20, phred = 33, 
    keepSortedBAM = FALSE, keepUniqBAM = FALSE)

snpi = snplists[matchsnp[i]]
shtid[i]
snpi
snpi = sprintf("%s/%s",snp.dir,snpi)


lngid = unlist(strsplit(fls,split="_filtered.bam"))
hap1 = sprintf("%s/%s_sorted_by_name_uniq_filtered_hap1.bam",outpath[i],lngid[i])
hap2 = sprintf("%s/%s_sorted_by_name_uniq_filtered_hap2.bam",outpath[i],lngid[i])
hapN = sprintf("%s/%s_sorted_by_name_uniq_filtered_hapN.bam",outpath[i],lngid[i])
if(!file.exists(hap1)){
  out = extractAsReads(input=louti, snpList=snpi, outputTag = NULL,prop.cut = 0.5) 
}else{
  out = c(hap1,hap2,hapN,louti)
}

#countReads(bamFile=louti, bedFile=, outputFile, overlapFraction = 0.9)

#do overall count for these files: hap1, hap2, hapN and total
counts=rep(0,4)
for(j in 1:4){
  system(sprintf("samtools view -c %s > tmp%s.txt",out[j],i))
  counti = read.table(sprintf("tmp%s.txt",i))
  counti
  counts[j] = counti[[1]]
}
counts

write.table(counts,sprintf("counts_%s_%s.txt",i,shtid[i]),row.names=F,quote=F,col.names=F)
file.remove(sprintf("tmp%s.txt",i))

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")

#HG38 anno
  txdb = loadDb(sprintf("%s/Homo_sapiens_GRCh38.sqlite", bow2.dir))
  seqlevels(txdb)
  columns(txdb)
  keytypes(txdb)
  
  genes = exonsBy(txdb, by="gene")

  bamfiles  = BamFileList(out, yieldSize=1000000)
  bamfiles
  

  se = summarizeOverlaps(features=genes, reads=bamfiles, mode="Union",
               singleEnd=FALSE, ignore.strand=FALSE, fragments=TRUE )
  
  se
  
  head(assay(se))
  colSums(assay(se))
  colData(se)
  rowData(se)
  str(metadata(rowData(se)))
  
  as1 = assay(se)

  
  write.table(as1, file = sprintf("%s/fl%s_%s_gene_level_counts.txt",count.dir,i,lngid[i]), append = FALSE,
    quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,
    col.names = TRUE)


q("no")
