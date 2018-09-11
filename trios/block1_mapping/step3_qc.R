args = commandArgs(TRUE)
args
k = as.numeric(args[1])
#k = 1
#sbatch --mem=16G --partition=day --wrap='R CMD BATCH step4_qc.R'
# Install R package ShorRead
# source("http://bioconductor.org/biocLite.R")
# biocLite("ShortRead", lib="~/R/Rlibs/")

#k = 1
library(ShortRead)

# following the workflow of Using Bioconductor for Sequence Data
# http://www.bioconductor.org/help/workflows/high-throughput-sequencing/


root.dir = "/lustre/scr/z/h/zhabotyn/R01"
pipe.dir = sprintf("%s/pipe", root.dir)
data.dir = sprintf("%s/data", root.dir)
data.dirj   = sprintf("%s/%s", data.dir, c("11_03_14", 
              "12_19_14", "12_29_14", "01_12_15", "01_22_15"))

for(k  in 1:length(data.dirj)){
#k = 5
  maxb = 151
  if(k==5)maxb = 76
  
  fls = list.files(data.dirj[k], "fastq$", full=TRUE)
  names(fls) = sub(".fastq", "", basename(fls))
  
  # generate qa summary
  date()
  qaSummary = qa(fls, type="fastq")
  date()
  
  dest.dir = sprintf("out%s/fast_qc_report",k)
  if(!file.exists(dest.dir))system(sprintf("mkdir %s",dest.dir))
  r12 = report(qaSummary, dest=dest.dir)
  
  rQ = qaSummary[["readQualityScore"]]
  dim(rQ)
  rQ[1:2,]
  
  perCycle = qaSummary[["perCycle"]]
  lapply(perCycle, dim)
  
  perCycle$baseCall[1:20,]
  perCycle$quality[1:20,]
  
  table(perCycle$quality$Cycle)
  table(perCycle$quality$lane)
  table(perCycle$quality$Score)
  
  sams = as.character(unique(perCycle$quality$lane))
  length(sams)
  
  mean = median = upper = lower = matrix(nrow=maxb, ncol=length(sams))
  colnames(mean) = sub(".fastq", "", sams)
  
  for(i in 1:length(sams)){
    
    sam1 = sams[i]
    dat1 = perCycle$quality[perCycle$quality$lane == sam1,]
    
    ind = 1:maxb
    if(k %in% c(3))ind = 1:76
    for(cycle in ind){
      dat2 = dat1[dat1$Cycle == cycle,]
      mean[cycle, i] = sum(dat2$Score*dat2$Count)/sum(dat2$Count)
      csum = cumsum(dat2$Count)
      
      nn  = csum[length(csum)]
      w25 = which(0.25*nn > csum[-length(csum)] & 0.25*nn <= csum[-1])
      w50 = which(0.50*nn > csum[-length(csum)] & 0.50*nn <= csum[-1])
      w75 = which(0.75*nn > csum[-length(csum)] & 0.75*nn <= csum[-1])
  
      median[cycle, i] = dat2$Score[w50+1]
      lower[cycle, i]  = dat2$Score[w25+1]
      upper[cycle, i]  = dat2$Score[w75+1]
  
    }
  }
  
  
  cols = rainbow(20)
  ltys = rep(1, 20)
  wwR2 = grep("R2", sams)
  ltys[wwR2] = 2
  
  summary(as.numeric(mean))
  pdf(sprintf("%s/perCycleScore_mean.pdf",dest.dir), width=8, height=4)
  par(mar=c(5,4,2,1), bty="n")
  plot(c(1,maxb), c(27,37), type="n", xlab="Cycle",
    ylab="Quality Score", main="mean")
  abline(v=seq(5,150,by=5), lty=2, col="grey")
  abline(h=seq(27,36,by=1), lty=2, col="grey")
  for(i in 1:length(sams)){
    lines(1:maxb, mean[,i], col=cols[i], lty=ltys[i])
  }
  legend("topright", c("R1", "R2"), lty=ltys, bg="white", horiz=TRUE)
  dev.off()
  
  sort(mean[50,])
  
  summary(as.numeric(median))
  pdf(sprintf("%s/perCycleScore_median.pdf",dest.dir), width=8, height=4)
  par(mar=c(5,4,2,1), bty="n")
  plot(c(1,maxb), c(31,38), type="n", xlab="Cycle",
  ylab="Quality Score", main="median")
  abline(v=seq(5,maxb-1,by=5), lty=2, col="grey")
  abline(h=seq(31,38,by=1), lty=2, col="grey")
  for(i in 1:length(sams)){
    lines(1:maxb, median[,i], col=cols[i], lty=ltys[i])
  }
  legend("bottom", c("R1", "R2"), lty=ltys, bg="white", horiz=TRUE)
  dev.off()
  
  summary(as.numeric(lower))
  pdf(sprintf("%s/perCycleScore_25.pdf",dest.dir), width=8, height=4)
  par(mar=c(5,4,2,1), bty="n")
  plot(c(1,maxb), c(12,37), type="n", xlab="Cycle",
  ylab="Quality Score", main="25 percentile")
  abline(v=seq(5,maxb-1,by=5), lty=2, col="grey")
  abline(h=seq(12,36,by=1), lty=2, col="grey")
  for(i in 1:length(sams)){
    lines(1:maxb, lower[,i], col=cols[i], lty=ltys[i])
  }
  legend("bottom", c("R1", "R2"), lty=ltys, bg="white", horiz=TRUE)
  dev.off()
  message("finished",k)
}


q("no")


