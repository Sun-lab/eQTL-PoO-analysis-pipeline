args = commandArgs(TRUE)
#args = c("2")
j = as.numeric(args[1])

get_block = function(str,split="_",block=1){
  unlist(strsplit(str,split=split))[block]
}
#ensure that when we write positions to a file R doesn't switch them to scientific
options("scipen"=999,"digits"=4)
#note that we also need to convert these snps to build 38
#I have converted by liftover positions of those snps

root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/omnidata", root.dir)
info.dir = sprintf("%s/info", root.dir)
pipe.dir = sprintf("%s/omnipipe", root.dir)
shape.dir = sprintf("%s/shapeit", root.dir)

geno.dir = sprintf("%s/data_genotype/omnireform", data.dir)
phas.dir = sprintf("%s/phasedR",geno.dir)
refinfo = sprintf("%s/1000GP_Phase3",info.dir)
impu.dir = sprintf("%s/imputed",geno.dir)

out.dir = sprintf("%s/snps_omni",data.dir)
if(!file.exists(out.dir))dir.create(out.dir)
lift.dir = sprintf("%s/liftover",root.dir)
vcf.dir = sprintf("%s/vcf",data.dir)
if(!file.exists(vcf.dir))dir.create(vcf.dir)

fls = list.files(impu.dir)
hap = fls[grep("hap",fls)]

alid = read.table(sprintf("%s/phasing_chr1.log.sample",phas.dir),as.is=T,skip=2)[,2]

smpls = read.table(sprintf("%s/prep/samples.dat",pipe.dir),as.is=T,header=T)[,1]
smpli = match(smpls,alid)
all(smpli==sort(smpli))
 
liftover.dir = sprintf("%s/liftover",root.dir)
setwd(liftover.dir)

#define snps (I'll exclude indels with matching to the base
bases=c("A","C","G","T")
maf = 0.01

#get all the appropriate chunks for the given chromosome
if(F){
fls = list.files(impu.dir)
to.rm = union(union(union(grep(".txt",fls),grep("info",fls)),grep("warnings",fls)),grep("summary",fls))
fls = fls[-to.rm]
hps = fls[grep("_haps",fls)]
alp = fls[grep("_allele_probs",fls)]
fls = setdiff(fls,union(hps,alp))
fls
hps
alp
}
#
off=5
#phe
#j=match("phased_imputed_chr2_85e6_90e6",fls)+1
#j = 1

   hpj = read.table(sprintf("%s/%s",impu.dir,hap[j]),as.is=T)
   hpj$V4 = as.character(hpj$V4)
   hpj$V5 = as.character(hpj$V5)
 
   rm.indels =  which(hpj$V4%in%bases & hpj$V5%in%bases)
   hpj = hpj[rm.indels,]
   dim(hpj)
   indloc = off +(smpli-1)*2+1
   indloc2 = off +(smpli-1)*2+2
   ind0 = sort(c(indloc,indloc2))
   ind = c(t(cbind(indloc,indloc2)))
   all(ind==ind0)
   inf = hpj[,1:5]   
   dat = hpj[,ind]
   to.kp = rowMeans(dat)
#   summary(to.kp[(to.kp>maf&to.kp<(1-maf))])
   to.kp = (to.kp>maf&to.kp<(1-maf))
   inf = inf[to.kp,]
   dat = dat[to.kp,]
   
   to.rm = as.numeric(get_block(hap[j],block=4))+5e6
   to.rm = which(inf[,3]==to.rm)
   if(length(to.rm>0)){
     inf = inf[-to.rm,]
     dat = dat[-to.rm,]
     message("removed boundary snp")
   }
   #
   chr = get_block(hap[j],block=3)
   nameit = which(inf[,2]==".")
   inf[nameit,2] = sprintf("%s:%s:%s:%s",chr,inf[nameit,3],inf[nameit,4],inf[nameit,5])
   lif = cbind(chr,start=inf[,3],end=inf[,3]+1,a1=inf[,4],a2=inf[,5],id=inf[,2])
   tmpj = sprintf("tmp_%s.bed",j)
   write.table(lif[,1:3],file=tmpj,sep=" ",
                quote=F,row.names=F,col.names=F)
    
   inp = sprintf("./liftOver %s/%s", liftover.dir,tmpj)
   cha = sprintf("%s/chains/hg19ToHg38.over.chain",liftover.dir)
   out = sprintf("%s/tmp_%s_hg38.bed",liftover.dir,j)
   unl = sprintf("%s/unm_%s_hg38.bed",liftover.dir,j)
   inp
   cha
   out
   unl
   comm = sprintf("%s %s %s %s",inp, cha,out,unl)
   comm
   system(comm)

   lif2 = read.table(out)
   info = file.info(unl)

   
   if(info$size>0){
     unmi2 = read.table(unl)
     to.rm = match(unmi2[,2],lif[,2])
     if(length(to.rm)>0)lif = lif[-to.rm,]
   }
file.remove(sprintf("%s/%s",liftover.dir,tmpj))
file.remove(out)
file.remove(unl)
   dim(lif)
   dim(lif2)
   lif[,2] = lif2[,2]    
   #once liftover is completed we again may have a situation with multiple snps at one position
   #if they don't conflict keep one of them, if they conflict, remove all
   if(nrow(lif)!=length(unique(lif[,2]))){
      to.fix = names(table(lif[,2])[table(lif[,2])>1])
      to.rm = numeric(0)
      for(fixi in to.fix){
        ambig = which(lif[,2]==fixi)
        attempt = aggregate(rep(0,length(ambig)),
        by=list(lif[ambig,1],lif[ambig,2],lif[ambig,5],lif[ambig,6]),FUN=sum)[,-5]
        if(nrow(attempt)==1){
          to.rm = c(to.rm,ambig[-1])
        }else{
          to.rm = c(to.rm,ambig)
        }
      }
      lif=lif[-to.rm,]
      message("removed ",length(to.rm)," ambiguous snps")
   }
   nrow(lif)==length(unique(lif[,2]))
   m = match(lif[,6],inf[,2])
   all(m==sort(m))
   inf = inf[m,]
   dat = dat[m,]

#chr21 9427957 chr21:9427957 T A . PASS AF=0.5506 GT      
   af = rowMeans(dat)
   lif = cbind(lif,af)
   lif2 = cbind(lif[,c(1,2,6,4,5)],".","PASS",round(af,4),"GT")
   smpli = 1:length(smpls)
   indloc = (smpli-1)*2+1
   indloc2 = (smpli-1)*2+2   
   vcf = matrix(".",nrow=nrow(dat),ncol=length(smpls))
   for(indi in smpli){
     vcf[,indi] = paste(dat[,indloc[indi]],dat[,indloc2[indi]],sep="|")
     #if(indi%%50==0)message(indi)
   }
   out = cbind(lif2,vcf)
vcfj = gsub("_haps","",hap[j])
vcfj = gsub("phased_imputed_","",vcfj)
write.table(out,sprintf("%s/SNP_VCF_%s.vcf",vcf.dir,vcfj),quote=F,row.names=F,col.names=F)

if(j==1){
  str = matrix(c("##fileformat=VCFv4.0",
  "##fileDate=20150416",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT "))
  str[4] = paste(str[4],paste(smpls,collapse = " "),sep = " ")
  write.table(str,sprintf("%s/templ.vcf",vcf.dir),quote=F,row.names=F,col.names=F)
}

q("no")
