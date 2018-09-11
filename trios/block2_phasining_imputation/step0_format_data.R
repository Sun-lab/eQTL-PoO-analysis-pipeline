args = commandArgs(TRUE)
args
chromosome = as.numeric(args[1])

#ensure that when we write positions to a file R doesn't switch them to scientific
options("scipen"=999,"digits"=4)

#for now I'll do build 37 (hg19) to match with 1000G which is in build 37
#when we create snps for individuals it would be needed to move them to build 38 (hg38) 
#in which the mapping of individuals was performed

#for PED file:
#do it for each chromosome separately:
#all our individuals are trios, so I'll get for each line a triple families, if 13281 is a family,
#(unless father is himself in a family)
#sex = 1 for male, =2 for female
#13281P NA12347 0 0 1 0 snps from genotype_auto_paternal.txt
#13281M NA12348 0 0 2 0 snps from genotype_auto_maternal.txt
#13281  NA12344 NA12347 NA12348 1 0 genotype_auto_individual.txt
#etc
#for MAP file:
#chr (int)
#snpid (str)
#snpgenpos (float) - put 0?
#snp physical pos(bp)

#chromosome = 10
root.dir  = "/home/groups/projects/Sun_RNA_seq"
info.dir  = sprintf("%s/info",root.dir  )
geno.dir  = sprintf("%s/data_genotype",root.dir)
lift.dir  = sprintf("%s/liftover",root.dir)

#get the list of how samples are related to each other
samples = read.table(sprintf("%s/samples2use.txt",info.dir),sep="\t",as.is=T,header=T)
#and  the sample ids from 1000 genomes
g1000 = read.table(sprintf("%s/1000IDS.txt",info.dir),as.is=T,sep="\t",header=T)

#check how many of samples of each kind are also in 1000G population
#Individual.ID,Paternal.ID,Maternal.ID:
table(!is.na(match(samples[,2],colnames(g1000))))
table(!is.na(match(samples[,3],colnames(g1000))))
table(!is.na(match(samples[,4],colnames(g1000))))

#do we have some sample to be a child of someone and a parent to someone else
any(samples[,2]%in%samples[,3])
any(samples[,2]%in%samples[,4])
any(samples[,3]%in%samples[,4])

indloc = sprintf("%s/genotype_auto_individuals.txt",geno.dir)
ind = read.table(indloc,sep="\t",header=T,as.is=T)
patloc = sprintf("%s/genotype_auto_paternal.txt",geno.dir)
pat = read.table(patloc,sep="\t",header=T,as.is=T)
matloc = sprintf("%s/genotype_auto_maternal.txt",geno.dir)
mat = read.table(matloc,sep="\t",header=T,as.is=T)
dim(ind)
dim(pat)
dim(mat)

#check whether we use same snps in all samples
#(if we had mismatches further step of taking intersection would be needed)
all(ind[,1]==pat[,1])
all(ind[,1]==mat[,1])

#check whether the order of parents in the data file matches
#the order of parents in family relationship file samples
all(match(samples[,2],colnames(ind))==match(samples[,3],colnames(pat)))
all(match(samples[,2],colnames(ind))==match(samples[,4],colnames(mat)))

#create location of the reformatted data (chromosome level PED/MAP)
reform.dir = sprintf("%s/ped90reform",geno.dir)
if(!file.exists(reform.dir))dir.create(reform.dir)


#load original hg18 positions file (or create it from the ind file
#this file can be directly used by liftover
hg18loc = sprintf("%s/snps_hg18.bed",info.dir)
hg19loc = sprintf("%s/snps_hg19.bed",info.dir)
un19loc = sprintf("%s/unm_hg19.bed",info.dir)
flag = file.exists(hg18loc)
if(flag){
  hg18 = read.table(hg18loc,as.is=T)
}else{
  #chr = unlist(strsplit(ind[,3],split="chr"))
  chr = as.character(ind[,3])
  chr2 = matrix(unlist(strsplit(ind[,3],split="chr")),nrow=2)[2,]
  starts = as.numeric(as.character(ind[,4]))
  ends = starts+1
  name = as.character(ind[,1])
  hg18 = data.frame(chr=chr,start=starts,end=ends,name=name)
  
  write.table(hg18,hg18loc,quote=F,sep=" ",row.names=F,col.names=F)
}

#converting to hg19
#note, liftOver should have proper permissions
#also, though I've seen some R related workflow
#www.bioconductor.org/help/workflows/liftOver
#it didn't install at the time, might change to it lateron 
flag = file.exists(hg19loc)
if(!flag){
  chain = sprintf("%s/chains/hg18ToHg19.over.chain.gz",lift.dir)
  #liftover format is just:
  #liftOver oldFile map.chain newFile unMapped
  lift = sprintf("%s/liftOver %s %s %s %s",
                 lift.dir,hg18loc,chain,hg19loc,un19loc)
  system(lift)
}

hg19 = read.table(sprintf(hg19loc,info.dir),as.is=T)
colnames(hg19) = colnames(hg18) = c("chr","start","end","name")

#if after lifting over two snps get to the same position
#we observer 3 such snps, I delete them both if they conflict, or duplicate
hg = apply(hg19[,1:2],1,paste,collapse=":")
hgt = table(hg)
table(hgt)
hgi = names(hgt)[hgt>1]
ind[!is.na(match(hg,hgi)),]
to.rm = numeric(0)
for(k in hgi){
  hgr = which(!is.na(match(hg,k)))
  cand = ind[hgr,]
  keep = TRUE  
  for(dbl in 2:nrow(cand)){
    keep = keep & all(unlist(cand[1,-c(1,2,4,6)])==unlist(cand[dbl,-c(1,2,4,6)]))
  }
  if(keep){
    to.rm = c(to.rm,hgr[-1])
  }else{
    to.rm = c(to.rm,hgr)
  }
}
to.rm
hg19 = hg19[-to.rm,]
  
dim(hg18)
dim(hg19)
nrow(hg18)-nrow(hg19)
hg18[1:5,]
hg19[1:5,]

#since we removed some snps we will have to remove corresponding data
keep = !is.na(match(hg18$name,hg19$name))
table(keep)
hg18 = hg18[keep,]
dim(hg18)

keep = !is.na(match(ind$rs,hg19$name))
ind = ind[keep,]
pat = pat[keep,]
mat = mat[keep,]
dim(ind)
dim(pat)
dim(mat)
c(all(hg19$name==ind$rs),all(hg19$name==pat$rs),all(hg19$name==mat$rs))

#process files into chromosome level
#note, if chromosome has being supplied script will process one chromosome
#otherwise it will loop over all autosomes

if(is.na(chromosome)){
  chrs = 1:22
}else{
  chrs = chromosome
}
for(i in chrs){
  chri = sprintf("chr%s",i)
  chri
  chri.dir = sprintf("%s/chri",reform.dir)
  if(!file.exists(chri.dir))dir.create(chri.dir)
  
  keep = hg19$chr==chri
  hg19i = hg19[keep,]
  indi = ind[keep,]
  pati = pat[keep,]
  mati = mat[keep,]
  dim(indi)
  dim(pati)
  dim(mati)

  #commenting this part, but may use in future  - removing monomorphic snps
  if(FALSE){
    indt = cbind(matrix(as.character(as.matrix(indi[,-(1:6)])),nrow=nrow(pati)),
                matrix(as.character(as.matrix(mati[,-(1:6)])),nrow=nrow(pati)),
                matrix(as.character(as.matrix(pati[,-(1:6)])),nrow=nrow(pati)))
    checkM = apply(indt,1,unique)
    check = function(list){
      sum(!list%in%"NN")>1
    }              
    keep = unlist(lapply(checkM,check))
    message("keep:",sum(keep)," remove:",sum(!keep))
    hg19i = hg19i[keep,]
    indi = indi[keep,]
    pati = pati[keep,]
    mati = mati[keep,]
  }
    
  ord = order(as.numeric(as.character(hg19i$start)))
  message("order in check: ",all(1:length(ord)==ord))
  hg19i = hg19i[ord,]
  indi = indi[ord,]
  pati = pati[ord,]
  mati = mati[ord,]

  #here goes PED file
  #in our data missing is N, in PED file it is expected to be 0
  prep = cbind(samples[,3],samples[,3],0,0,1,0)
  prem = cbind(samples[,4],samples[,4],0,0,2,0)
  prei = cbind(samples[,2],samples[,2],samples[,3],samples[,4],samples[,5],0)
  
  pato = unlist(strsplit(as.character(as.matrix((pati[,match(prep[,2],colnames(pati))]))),""))
  pato[pato=="N"] = "0";table(pato)
  pato = matrix(pato,nrow=2*nrow(pati))
  
  mato = unlist(strsplit(as.character(as.matrix((mati[,match(prem[,2],colnames(mati))]))),""))
  mato[mato=="N"] = "0";table(mato)
  mato = matrix(mato,nrow=2*nrow(mati))
  
  indo = unlist(strsplit(as.character(as.matrix((indi[,match(prei[,2],colnames(indi))]))),""))
  indo[indo=="N"] = "0";table(indo)
  indo = matrix(indo,nrow=2*nrow(indi))
  
  pat0 = cbind(prep,t(pato))
  mat0 = cbind(prem,t(mato))
  ind0 = cbind(prei,t(indo))
  all0 = t(matrix(c(t(cbind(pat0,mat0,ind0))),nrow=ncol(pat0)))
  all0[1:12,1:15]
  for(m in 1:nrow(pat0)){
    message(all(all0[3*(m-1)+1,]==pat0[m,])," ",all(all0[3*(m-1)+2,]==mat0[m,]),
           " ",all(all0[3*(m-1)+3,]==ind0[m,]))
  }
  write.table(all0,sprintf("%s/geno_chr%s.ped",chri.dir,chrs[i]),
                   row.names=F,col.names=F,quote=F)
  
  #this is MAP file
  hg19i$chr = i
  #I reuse end column to define absent phenotype, that is a required field
  hg19i$end = 0
  hg19i = hg19i[,c(1,4,3,2)]
  write.table(hg19i,sprintf("%s/geno_chr%s.map",chri.dir,chrs[i]),
                    row.names=F,col.names=F,quote=F)
  message(i)
#}

q("no")
