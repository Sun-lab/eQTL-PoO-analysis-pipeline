#imputation
#ensure that when we write positions to a file R doesn't switch them to scientific
options("scipen"=999,"digits"=4)

queue = "bat" 

root.dir = "/lustre/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/data", root.dir)
info.dir = sprintf("%s/info", root.dir)
pipe.dir = sprintf("%s/pipe", root.dir)
shape.dir = sprintf("%s/shapeit", root.dir)

geno.dir = sprintf("%s/data_genotype/ped90reform", data.dir)
phas.dir = sprintf("%s/phasedR",geno.dir)
refinfo = sprintf("%s/1000GP_Phase3",info.dir)
if(!file.exists(phas.dir))dir.create(phas.dir)
refinfo = sprintf("%s/1000GP_Phase3",info.dir)

impute.dir = sprintf("%s/impute2/", root.dir)
impu.dir = sprintf("%s/imputed",geno.dir)
#setwd(impute.dir)
#or
#module add impute if available

refinfo = sprintf("%s/1000GP_Phase3",info.dir)
if(!file.exists(impu.dir))dir.create(out.dir)

#i = 22;j = 2
#this file will keep the split positions for each chromosome just in case they will be needed later
write.table(data.frame("chr","ind"),sprintf("%s/step5_indicies.txt",pipedir),
            append=F,quote=F,row.names=F,col.names=F,sep="\t")
#we still use the same 1000G reference:
all1000G = "1000GP_Phase3" #common reference prefix
mem = "32G"
#processing  each chromosome, impute runs on 5MB chunks (no more than 7MB), so first split into chunks
for(i in 1:22){
  #using the appropriate positions from both data and 1000G we find the range in which we will do imputation
  khf = sprintf("%s/phasing_chr%s.log.haps",phas.dir,i)
  dati = read.table(khf,as.is=T)
  rng = range(dati[,3])
  
  chri = sprintf("chr%s",i)
  con = gzfile(sprintf("%s/%_%s.legend.gz",refinfo,all1000G,chri))
  refi = readLines(con)
  close(con)
  rng2 = range(as.numeric(matrix(unlist(strsplit(refi[-1]," ")),nrow=12)[2,]))
  rng[1] = min(rng[1],rng2[1])
  rng[2] = max(rng[2],rng2[2])
  message(paste(rng2-rng,collapse=" "))
  
  mb = 1e6
  indj = sprintf("%se6",c(seq(floor(rng[1]/1e6),floor(rng2[2]/mb),by=5),ceiling(rng[2]/mb)))
  write.table(data.frame(i,paste(indj,collapse=";")),sprintf("%s/indicies.txt",pipe.dir),
              append=T,quote=F,row.names=F,col.names=F,sep="\t")
  #since we prephased using shapeit we need to use -use_prephased_g and -known_haps_g
  pref = "./impute2 -use_prephased_g"
  m = sprintf("-m %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)
  h = sprintf("-h %s/%s_chr%s.hap.gz",refinfo,all1000G,i)
  l = sprintf("-l %s/%s_chr%s.legend.gz ",refinfo,all1000G,i)
  khg = sprintf("-known_haps_g %s",khf)
  oth = "-align_by_maf_g -Ne 20000 -seed 12345"

  #here we will run by the 5MB increment from the rng[1] until we exaust the chromosome
  for(j in 2:length(indj)){
    #j = length(indj)
    out = sprintf("-phase -o %s/phased_imputed_chr%s_%s_%s",impu.dir,i,indj[j-1],indj[j])    
    pos = sprintf("-int %s %s",indj[j-1]+1,indj[j])

    com = sprintf("%s %s %s %s %s %s %s %s",pref,m,h,l,khg,oth,pos,out)
  qout = sprintf("--output=%s/rout/step5_out_%s.out",pipe.dir,i)    
  com2 = sprintf("sbatch %s --partition=%s --mem=%s --wrap='%s'",qout,queue,mem,com)
    message(com2)
    system(com2)
  }
}

q("no")
