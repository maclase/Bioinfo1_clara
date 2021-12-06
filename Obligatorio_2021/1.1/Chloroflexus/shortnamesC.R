filename="ChloroflexusCDS" # el archivo original, en este caso, se llama Ecoli.cds

library(seqinr)
fasnuc = read.fasta(paste(filename,"pep",sep = ".")) # el nombre real del archivo
FASNUC <- sapply(fasnuc, toupper)
len = length(FASNUC)
head(names(FASNUC))
namesOrig = names(FASNUC)

namesNew=c()
for(i in 1:len){
  namesNew[i]= paste("c",i,sep=".")
}

head(namesNew)



write.fasta(FASNUC,names=namesNew,file.out = paste("nc",filename,"PEP",sep=".")) 
# nc =nombres cortos
faspep = lapply(FASNUC,translate)
write.fasta(faspep,names = namesNew,file.out = paste("nc",filename,"pep",sep="."))

dfnames=data.frame(namesNew,namesOrig)
head(dfnames)
write.table(dfnames,quote=FALSE,row.names = FALSE,
            file=paste("nc",filename,sep = "."),sep="\t")


###### fin del script (lo que sigue es para hacer cosas en clase) #######



system("./runBlastp.bash nc.NeisseriaCDS.pep nc.argN.PEP")

t=read.table("tablast.nc.NeisseriaCDS.pep-nc.argN.PEP")
colnames(t)= c("query","subject","score","len_aln","mism","gapop",
      "qstart","qend","sstart","send","evalue","score","qlen","slen")
head(t)

# Podemos filtrar la tabla desde aquí

order(t)

borrar=which(t$query == t$subject)
?which
t2=t[-borrar,]
dim(t)
dim(t2)

borrar=which(t2$p_id_aln<50)
t3 = t2[-borrar,]
dim(t3)

#makeblastdb -in nc.argN.PEP -dbtype prot -out argN_db
#blastp -query nc.NeisseriaCDS.pep -db argN_db -out tablast_Neiss_argN -outfmt "6 std qseqid sseqid qlen slen length pident evalue bitscore"