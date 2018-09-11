#for one file run tophat with option -G

#module add bowtie2
#chmod 755 *
root.dir = "/lustre/scr/z/h/zhabotyn/R01"
pipe.dir = sprintf("%s/omnipipe", root.dir)
data.dir = sprintf("%s/omnidata", root.dir)

gtf.dir = sprintf("%s/bowtie2_indexes", root.dir)
gtf.file = sprintf("%s/gencode.v21.annotation.gtf", gtf.dir)
tra =sprintf("--transcriptome-index=%s/genome/gencode.v21.annotation", gtf.dir)
idx = sprintf("%s/genome/hg38", gtf.dir)


subi = 1

  fastq.dir1 = file.path(data.dir)
  fastq.dir2 = file.path(data.dir, "fastq_trimmed_and_filtered")
  output.dir = file.path(data.dir, "tophat_output_liberal")
  if(!file.exists(output.dir))system(sprintf("mkdir %s",output.dir))
  output.file = sprintf("step6_tophat_liberal_0.sh")
  batch=sprintf("%s/R_batch",pipe.dir)
  if(!file.exists(batch))dir.create(batch)
  
  
  setwd(fastq.dir1)
  
  fls = list.files(fastq.dir1, "fastq.gz$", full=TRUE)
  names(fls) = sub(".fastq.gz", "", basename(fls))
  samples    = gsub("_1", "", names(fls))
  samples    = unique(gsub("_2", "", samples))
  samples
  
  setwd(batch)
  
  cat("", file = output.file)
  
#for one file run with option -G
  
  ind = 1:length(samples)

    samp = samples[subi]

    cmd0 = "sbatch -n 8  --partition=week --wrap=\'tophat -G"    
    cmd1 = sprintf("%s %s %s -o %s/%s", cmd0, gtf.file, tra, output.dir, samp)    

    cmd1 = sprintf("%s --library-type=fr-firststrand --read-mismatches 3", cmd1)
    cmd1 = sprintf("%s --read-gap-length 3 --read-edit-dist 3", cmd1)
    cmd1 = sprintf("%s --read-realign-edit-dist 0 -p 8 %s", cmd1, idx)
    cmd1 = sprintf("%s %s/%s_R1_trimmed_filtered.fastq.gz", cmd1, fastq.dir2, samp)
    cmd1 = sprintf("%s %s/%s_R2_trimmed_filtered.fastq.gz", cmd1, fastq.dir2, samp)
    cmd1 = sprintf("%s,%s/%s_Si_trimmed_filtered.fastq.gz\'", cmd1, fastq.dir2, samp)
    
    system(cmd1)

q("no")

