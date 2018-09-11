library(Rsamtools)

root.dir = "/lustre/scr/z/h/zhabotyn/R01"
#root.dir = "G:/projects/Sun_RNA_seq"
data.dir = sprintf("%s/data", root.dir)
geno.dir = sprintf("%s/data_genotype", root.dir)
cnt.dir = sprintf("%s/counts/ind_new", data.dir)
#cnt.dir = sprintf("%s/counts/ind", data.dir)
tre.dir = sprintf("%s/counts/ind_conv",data.dir)
if(!file.exists(tre.dir))dir.create(tre.dir)
info.dir = sprintf("%s/info", root.dir)
setwd(tre.dir)

ind = read.table(sprintf("%s/genotype_X_individuals.txt", geno.dir), nrow=1)

get_block = function(x, split="-", block=2){
  unlist(strsplit(x, split=split))[block]
}
samples = list.files(cnt.dir, pattern="\\.txt")
pos = as.numeric(gsub("fl", "", sapply(samples, get_block, split="_", block=1)))
ids = sapply(sapply(samples, get_block), get_block, split="_", block=1)
ids = gsub("Gm","NA", ids)
o = order(pos)
bat = rep(1, 50); bat[pos>10] = 2; bat[pos>30] = 3;names(bat)=ids
ids = ids[o]
samples = samples[o]
uids = unique(ids)

for(i in 1:length(samples)){
  if(i == 1){
    cnts = read.table(sprintf("%s/%s", cnt.dir,samples[i]),sep="\t")
    tots = hap2 = hap1 = matrix(0, nrow=nrow(cnts), ncol=30)
    colnames(tots) = colnames(hap2) = colnames(hap1) = uids
    rownames(tots) = rownames(hap2) = rownames(hap1) = rownames(cnts)
  }else{
    cnti = read.table(sprintf("%s/%s",cnt.dir,samples[i]),sep="\t")
    if(!all(rownames(cnti)==rownames(cnts))){
      message("not matching genenames")
      break
    }
    cnts = cnti
  }
  m = match(ids[i], uids)
  hap1[,m] = hap1[,m] + cnts[,1]
  hap2[,m] = hap2[,m] + cnts[,2]
  tots[,m] = tots[,m] + cnts[,4]

  if(i%%5==0)message("processed ",i,"'th sample")
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
  kp = which(apply(hapi>4, 1, sum)>1)
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
#get log(trc), pc1, pc2, pc3, pc4, pc5, mon1, mon2 design matrix
#write to
#samples.dat



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
pchs = rep(1:3, each=10)
cols = rep(c("black", "blue", "gold"), each=10)
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
samdat$PC5 =  pcaA$x[,5]
samdat$gr1 = c(rep(0, 10), rep(1, 10), rep(0, 10))
samdat$gr2 = c(rep(0, 10), rep(1, 10), rep(0, 10))
write.table(samdat,sprintf("%s/samples.dat",tre.dir),row.names=F,col.names=T,sep="\t",quote=F)


q("no")
