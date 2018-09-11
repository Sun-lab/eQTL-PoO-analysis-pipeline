root.dir = "/lustre/scr/z/h/zhabotyn/R01"
#root.dir = "G:/projects/Sun_RNA_seq"
data.dir = sprintf("%s/omnidata", root.dir)
par.dir = sprintf("%s/counts", data.dir)
cnt.dir = sprintf("%s/cntu38", par.dir)
tre.dir = sprintf("%s/parplus_conv",data.dir)
if(!file.exists(tre.dir))dir.create(tre.dir)
info.dir = sprintf("%s/info", root.dir)
setwd(tre.dir)

samp = read.table(sprintf("%s/geuvadis_select.txt", par.dir), header=T)

#note, for this larger number of samples we choose genes having at least 5 samples with at least 5 allele specific counts

get_block = function(x, split="-", block=2){
  unlist(strsplit(x, split=split))[block]
}
samples = list.files(cnt.dir, pattern="\\.txt")
ids = sapply(samples, get_block, split="_")
m = match(samp[,1], ids)
samples = samples[m]
ids = ids[m]

for(i in 1:length(samples)){
  if(i == 1){
    cnts = read.table(sprintf("%s/%s", cnt.dir,samples[i]),sep="\t")
    tots = hap2 = hap1 = matrix(0, nrow=nrow(cnts), ncol=length(ids))
    colnames(tots) = colnames(hap2) = colnames(hap1) = ids
    rownames(tots) = rownames(hap2) = rownames(hap1) = rownames(cnts)
  }else{
    cnti = read.table(sprintf("%s/%s",cnt.dir,samples[i]),sep="\t")
    if(!all(rownames(cnti)==rownames(cnts))){
      message("not matching genenames")
      break
    }
    cnts = cnti
  }
  m = match(ids[i], ids)
  hap1[,m] = hap1[,m] + cnts[,1]
  hap2[,m] = hap2[,m] + cnts[,2]
  tots[,m] = tots[,m] + cnts[,4]

  if(i%%25==0)message("processed ",i,"'th sample")
}
geneinfo = read.csv(sprintf("%s/genepos_GRCh38.csv", info.dir), as.is=T)
o = order(geneinfo$chr, geneinfo$start)
geneinfo = geneinfo[o,]
m = match(geneinfo$id, rownames(cnts))
table(geneinfo$id==rownames(cnts)[m])
hap1 = hap1[m,]
hap2 = hap2[m,]
tots = tots[m,]

#get i'th chromosome and print it
i = 22
for(i in 1:22){
  chri = sprintf("chr%s", i)
  kp = geneinfo$chr==chri
  infoi = geneinfo[kp,]
  hap1i = hap1[kp,]
  hap2i = hap2[kp,]
  totsi = tots[kp,]
  hapi = hap1i + hap2i
  kp = which(apply(hapi>4, 1, sum)>4)
  infoi = infoi[kp,]
  hap1i = hap1i[kp,]
  hap2i = hap2i[kp,]
  totsi = totsi[kp,]

  write.table(hap1i, sprintf("gene_a1_chr%s.dat", i), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(hap2i, sprintf("gene_a2_chr%s.dat", i), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(totsi, sprintf("gene_Tcount_chr%s.dat", i), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(infoi, sprintf("gene_info_chr%s.dat", i), row.names=F, col.names=T, quote=F, sep="\t")
}
dim(hap1i)

#geneinfo
#id	chr	start	end
#finish this - write to gene_info_chr22.dat
#get log(trc), pc1, pc2, pc3, pc4
#write to samples.dat



#get pci
trc = tots
tct = colSums(trc)
out3 = log10(tct)
trc = trc%*%diag(median(tct)/tct)
dim(trc)
summary(colSums(trc))
trc = trc[apply(trc,1,max)>50,]
dim(trc)
summary(colSums(trc))

#pca
x = rep(1:3,each=10)
trc = log(trc+1)
pcaA = prcomp(t(trc),retx=T)
summary(pcaA)


if(!file.exists("figures"))dir.create("figures")
cx = 1
pchs = rep(1, nrow(samp))
pchs[samp$GBR==1] = 2
pchs[samp$FIN==1] = 3
cols = rep("black", nrow(samp))
cols[samp$GBR==1] = 2
cols[samp$FIN==1] = 3
bitmap(paste("figures/total_counts_PCA_50_alt.png",sep=""), width=6, height=6,res=300,units="in")
par(mar=c(5,5,1,0), mfrow=c(2,2))
  barplot(summary(pcaA)$importance[2,1:20],
  xlab="principal component",ylab="proportion",cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,at=1:20,labels=sprintf("PC%s",1:20),cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  
  tti = pcaA$x
  pc1i = 1;  pc2i = 2
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 1;  pc2i = 3
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC3",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 2;  pc2i = 3
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC2",ylab="PC3",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
dev.off()

bitmap(paste("figures/total_counts_PCA_50_4.png",sep=""), width=6, height=6,res=300,units="in")
par(mar=c(5,5,1,0), mfrow=c(2,2))
  pc1i = 1;  pc2i = 4
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 2;  pc2i = 4
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 3;  pc2i = 4
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
dev.off()

bitmap(paste("figures/total_counts_PCA_50_5.png",sep=""), width=6, height=6,res=300,units="in")
par(mar=c(5,5,1,0), mfrow=c(2,2))
  pc1i = 1;  pc2i = 5
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 2;  pc2i = 5
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 3;  pc2i = 5
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 4;  pc2i = 5
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
dev.off()

bitmap(paste("figures/total_counts_PCA_50_6.png",sep=""), width=6, height=9,res=300,units="in")
par(mar=c(5,5,1,1), mfrow=c(2,3))
  pc1i = 1;  pc2i = 6
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 2;  pc2i = 6
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 3;  pc2i = 6
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 4;  pc2i = 6
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
  pc1i = 5;  pc2i = 6
  plot(tti[,pc1i], tti[,pc2i],  bty="n", pch=pchs,col=cols,xlab="PC1",ylab="PC2",cex=cx,cex.lab=cx,xaxt="n",yaxt="n")
  axis(1,cex.lab=cx,cex.axis=cx)
  axis(2,cex.lab=cx,cex.axis=cx)
dev.off()


samdat = data.frame(ID=names(out3), ttlCount=out3)
samdat$PC1 =  pcaA$x[,1]
samdat$PC2 =  pcaA$x[,2]
samdat$PC3 =  pcaA$x[,3]
samdat$PC4 =  pcaA$x[,4]
samdat$GBR = samp$GBR
samdat$FIN = samp$FIN
samdat$TSI = samp$TSI
#note, we do see that PC1, PC2 and PC3 reflect whether people are GBR, FIN or TSI so in analysis use actual GBR, FIN and TSI grouping
write.table(samdat,sprintf("%s/samples.dat",tre.dir),row.names=F,col.names=T,sep="\t",quote=F)


q("no")
