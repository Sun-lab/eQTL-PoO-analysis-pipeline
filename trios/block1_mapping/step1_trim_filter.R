args = commandArgs(trailingOnly=TRUE)
args
subk = as.numeric(args[1])
subk
#subk = 1
#sbatch --mem=16G --partition=day --wrap='R CMD BATCH "--args 1" step5_trim_filter.R'
library(ShortRead)

myFilterAndTrim = function(samp, fastq.dir, output.dir)
{
  
  fl1 = sprintf("%s/%s_1.fastq.gz", fastq.dir, samp)
  fl2 = sprintf("%s/%s_2.fastq.gz", fastq.dir, samp)

  flo1 = sprintf("%s/%s_1_trimmed_filtered.fastq.gz", output.dir, samp)
  flo2 = sprintf("%s/%s_2_trimmed_filtered.fastq.gz", output.dir, samp)
  flo3 = sprintf("%s/%s_Si_trimmed_filtered.fastq.gz", output.dir, samp)
  if(file.exists(flo1))file.remove(flo1)
  if(file.exists(flo2))file.remove(flo2)
  if(file.exists(flo3))file.remove(flo3)

  ## open input stream
  stream1 = open(FastqStreamer(fl1))
  stream2 = open(FastqStreamer(fl2))

  j = 0

  repeat {
    j = j + 1
    
    ## input chunk
    fq1 = yield(stream1)
    fq2 = yield(stream2)

    # encoding(quality(fq1))
    if (length(fq1) != length(fq2)){
      stop("mismatch length beween two files\n")
    }
    
    if (length(fq1) == 0){
      break
    }
    
    id1 = as.character(id(fq1))
    id2 = as.character(id(fq2))
    
    id1m = matrix(unlist(strsplit(id1, split=" ")), ncol=2, byrow=TRUE)
    id2m = matrix(unlist(strsplit(id2, split=" ")), ncol=2, byrow=TRUE)

    dim(id1m)
    dim(id2m)
    id1m[1:2,]
    id2m[1:2,]

    if (any(id1m[,1] != id2m[,1])){
      stop("mismatch of sequence ids\n")
    }


    str1 = sprintf("%d: length(fq1)/length(fq2)=%d/%d\n",
                    j, length(fq1), length(fq2))
    cat(str1)
    
    ## trim and filter, e.g., reads cannot contain 'N'...
    flag1 = flag2 = TRUE
    fq1 = fq1[nFilter()(fq1)]
    fq2 = fq2[nFilter()(fq2)]
    if(length(fq1)==0){
      flag1=FALSE
      message("long set of reads with N's at set1: ",j)
    }
    if(length(fq2)==0){
      flag2=FALSE
      message("long set of reads with N's at set2: ",j)
    }
    
    str1 = "After removing reads containing 'N', length(fq1)/length(fq2)"
    str1 = sprintf("%s=%d/%d\n", str1, length(fq1), length(fq2))
    cat(str1)

    if(flag1){
    ## trim as soon as 4 of 5 nucleotides has quality encoding less
    ## than or equalt to "4" (phred score 19)
    fq1 = trimTailw(fq1, 4, "4", 2)
    ## drop reads that are less than 75bp
    fq1 = fq1[width(fq1) >= 75]
    }
    if(flag2){
    ## trim as soon as 4 of 5 nucleotides has quality encoding less
    ## than or equalt to "4" (phred score 19)
    fq2 = trimTailw(fq2, 4, "4", 2)
    ## drop reads that are less than 75bp
    fq2 = fq2[width(fq2) >= 75]
    }
    str1 = "After trimming & removing reads < 75bp, length(fq1)/length(fq2)"
    str1 = sprintf("%s=%d/%d\n\n", str1, length(fq1), length(fq2))
    cat(str1)
    
    if(length(fq1)==0){
      flag1=FALSE
      message("all remaining reads were removed in set1: ",j)
    }
    if(length(fq2)==0){
      flag2=FALSE
      message("all remaining reads were removed in set2: ",j)
    }
    
    id1 = as.character(id(fq1))
    id2 = as.character(id(fq2))
    
    id1m = sapply(strsplit(id1, split=" "), function(x){x[[1]]})
    id2m = sapply(strsplit(id2, split=" "), function(x){x[[1]]})

    iduse = intersect(id1m, id2m)

    mm1 = match(iduse, id1m)
    mm2 = match(iduse, id2m)
    
    fq1p = fq1[mm1]
    fq2p = fq2[mm2]

    if(length(mm1)>0){
      fq1s = fq1[-mm1]
    }else{
      fq1s = fq1
    }
    if(length(mm2)>0){
      fq2s = fq2[-mm2]
    }else{
      fq2s = fq2
    }

    ## append to destination
    writeFastq(fq1p, flo1, "a")
    writeFastq(fq2p, flo2, "a")
    writeFastq(fq1s, flo3, "a")
    writeFastq(fq2s, flo3, "a")
  }
  
  close(stream1)
  close(stream2)

}


root.dir = "/lustre/scr/z/h/zhabotyn/R01"
setwd(root.dir)
nrd = 1e6
data.dir   = sprintf("%s/data", root.dir)
data.dirj   = sprintf("%s/%s", data.dir, c("11_03_14", 
              "12_19_14", "12_29_14", "01_12_15", "01_22_15"))

  output.dir = file.path(data.dirj[subk], "fastq_trimmed_and_filtered")
  if(!file.exists(output.dir))dir.create(output.dir)

  setwd(data.dirj[subk])

  fls = list.files(data.dirj[subk], "fastq.gz$", full=TRUE)
  names(fls) = sub(".fastq.gz", "", basename(fls))
  samples    = gsub("_1", "", names(fls))
  samples    = unique(gsub("_2", "", samples))
  samples

  for(i in 1:length(samples)){
    cat("\n------------------------------------------------------------\n")
    cat(i, fls[i], date())
    cat("\n------------------------------------------------------------\n")
  
    samp = samples[i]
    myFilterAndTrim(samp, data.dirj[subk], output.dir)
  }

warnings()

q("no")

