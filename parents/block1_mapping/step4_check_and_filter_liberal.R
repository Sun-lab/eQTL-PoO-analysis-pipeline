args = commandArgs(TRUE)
args
subi = as.numeric(args[1])

library(Rsamtools)

# -----------------------------------------------------------------
# read in sample information
# -----------------------------------------------------------------

root.dir = "/lustre/scr/z/h/zhabotyn/R01"
pipe.dir = sprintf("%s/omnipipe", root.dir)
data.dir = sprintf("%s/omnidata", root.dir)

setwd(pipe.dir)

  dati = sprintf("%s/tophat_output_liberal/",data.dirj[k])
  bamFiles = list.files(path=dati,
    pattern="accepted_hits.bam$", recursive=TRUE)
  to.rm = grep("tmp",bamFiles)
  if(length(to.rm)>0)bamFiles = bamFiles[-to.rm]
  bamFiles
  
  unMapped = list.files(dati,
    pattern="unmapped.bam$",
  recursive=TRUE)
  to.rm = grep("tmp",unMapped)
  if(length(to.rm)>0)unMapped = unMapped[-to.rm]
  unMapped
    
  
  ctMtrix = NULL
  
  outi =  sprintf("%s/bam_liberal",dati)    
    i = subi

    cat(i, date(), "\n")
    
    ctv  = rep(0, 3)
    
    bami = file.path(dati, bamFiles[i])
    umpi = file.path(dati, unMapped[i])
  
    # ----------------------------------------------------------
    # index the bam file
    # ----------------------------------------------------------
  
    bamIdx = paste(bami, ".bai", sep="")
    if(!file.exists(bamIdx)){
      bamIdx = indexBam(bami)
    }
    
    umpIdx = paste(umpi, ".bai", sep="")
    if(!file.exists(umpIdx)){
      umpIdx = indexBam(umpi)
    }
  
    # ----------------------------------------------------------
    # counting
    # ----------------------------------------------------------
    
    date()
    ct1 = countBam(bami, index=bamIdx)
    date()
    ct2 = countBam(umpi, index=umpIdx)
    date()
  
    ctv[1] = ct1$records + ct2$records
    ctv[2] = ct1$records
  
    cat(round(ct1$nucleotides/ct1$records,2), " ")
    cat(round(ct2$nucleotides/ct2$records,2), " ")
  
    # ----------------------------------------------------------
    # getUnique and filtering
    # ----------------------------------------------------------
    
    flag1  = scanBamFlag(isUnmappedQuery=FALSE, isNotPrimaryRead=FALSE,
                        isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
    param1 = ScanBamParam(flag=flag1, what="seq")

    if(!file.exists(outi))system(sprintf("mkdir %s",outi))
    sampi  = sub(dati, "", dirname(bami))
    bamNi  = paste(sampi, "_filtered.bam", sep="")
    bamNi  = sprintf("%s%s",outi, bamNi)

    date()
    if(!file.exists(bamNi))filterBam(bami, destination=bamNi, index=bamIdx, param=param1)
    date()
  
    # ----------------------------------------------------------
    # counting again
    # ----------------------------------------------------------
    
    ct3 = countBam(bamNi, index=paste(bamNi, ".bai", sep=""))
    ctv[3] = ct3$records
    cat(round(ct3$nucleotides/ct3$records,2), " \n")
  
    ctMtrix = rbind(ctMtrix, c(sampi, ctv))
  }
  
  colnames(ctMtrix) = c("sample", "all", "mapped", "filtered")
  
  write.table(ctMtrix, file = sprintf("%s/counts_liberal_%s.txt",outi, subi),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


q(save="no")

