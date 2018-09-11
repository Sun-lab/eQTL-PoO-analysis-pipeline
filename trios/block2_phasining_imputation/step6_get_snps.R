args = commandArgs(TRUE)
#args = c("1")
chri = args[1]

#ensure that when we write positions to a file R doesn't switch them to scientific
options("scipen"=999,"digits"=4)
#note that we also need to convert these snps to build 38
#I have converted by liftover positions of those snps


root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/data", root.dir)
info.dir = sprintf("%s/info", root.dir)
pipe.dir = sprintf("%s/pipe", root.dir)
shape.dir = sprintf("%s/shapeit", root.dir)

geno.dir = sprintf("%s/data_genotype/ped90reform", data.dir)
phas.dir = sprintf("%s/phasedR",geno.dir)
if(!file.exists(phas.dir))dir.create(phas.dir)
refinfo = sprintf("%s/1000GP_Phase3",info.dir)
impu.dir = sprintf("%s/imputed",geno.dir)

out.dir = sprintf("%s/snps_ext",data.dir)
if(!file.exists(out.dir))dir.create(out.dir)
setwd(impu.dir)

#the list of samples for which we will get snps
samples = read.table(sprintf("%s/phasing_chr22.log.sample",phas.dir,chri),header=T,as.is=T)
samples = samples[-1,1]
#define snps (I'll exclude indels with matching to the base
bases=c("A","C","G","T")

#get all the appropriate chunks for the given chromosome
fls = list.files(pattern=sprintf("chr%s_",chri))
to.rm = union(union(union(grep(".txt",fls),grep("info",fls)),grep("warnings",fls)),grep("summary",fls))
fls = fls[-to.rm]
hps = fls[grep("_haps",fls)]
alp = fls[grep("_allele_probs",fls)]
fls = setdiff(fls,union(hps,alp))
fls
hps
alp

append=F
#
j = 1;off=5
probs = c(.9,.95,.99,.995,.999)

summ_het = matrix(0,nrow=90,ncol=length(probs))
rownames(summ_het) = samples
colnames(summ_het)=sprintf("p%s",probs)
tot = rep(0,90)

#phe
for(j in 1:length(fls)){   
#j=match("phased_imputed_chr2_85e6_90e6",fls)+1
#j = 1

   flj = fls[j]
   snpj = read.table(flj,as.is=T)
   hpj = read.table(hps[j],as.is=T)
   alj = read.table(alp[j],as.is=T)
   snpj$V4 = as.character(snpj$V4)
   snpj$V5 = as.character(snpj$V5)
 
   rm.indels =  which(snpj$V4%in%bases & snpj$V5%in%bases)
   snpj = snpj[rm.indels,]
   alj = alj[rm.indels,]
   hpj = hpj[rm.indels,]
   dim(hpj)
   dim(alj)
   dim(snpj)

   for(ind in 1:90){
   #ind = 1
   #ind=ind+1
     nmind=sprintf("ind%s_%s",ind,samples[ind])
     out = sprintf("%s/%s",out.dir,nmind)
     if(!file.exists(out))dir.create(out)
     hetloc = off+(ind-1)*3+2
     indloc = off +(ind-1)*2+1
     for(probk in 1:length(probs)){
       summ_het[ind,probk] = summ_het[ind,probk] + sum(snpj[,hetloc] >= probs[probk])
     }
     tot[ind] = tot[ind] + nrow(snpj)
     #I'll use only the snps with probability>0.999
     flag = snpj[,hetloc]>probs[5]
     if(sum(flag)>0){
       snpind = snpj[flag,]
       alind = alj[flag,c(1:5,indloc,indloc+1)]
       #note, since we select only heterozygous they are either in the same order as in columns 4,5 
       #(if code 0,1 for appropriate individual) or a flip of columns 4,5 (if code is 1,0)
       hpind0 = hpind = hpj[flag,c(1:5,indloc,indloc+1)]
       flip = (hpind[,6]==1)
       isalt = as.numeric(flip)
       hpind[,1] = sprintf("chr%s",chri)
       hpind[flip,4:5] = hpind[flip,5:4]
       hpind = cbind(hpind, isalt)
       outfile = sprintf("%s/chr%s.txt",out,chri)
       app=file.exists(outfile)
       write.table(hpind[,c(1,3,4,5,2,8)],outfile,row.names=F,col.names=F,quote=F,append=app)
     }
   }

   message(flj)
}

#summary of 
outfile = sprintf("%s/chr_%s_num_snps_1.txt",out.dir,chri)
write.table(cbind(summ_het,tot),outfile,row.names=T,col.names=T,quote=F)

q("no")
