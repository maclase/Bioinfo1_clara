library(ade4)
s=read.fasta("Pmarinus.ffn")
length(s)
ucodones=sapply(s,uco)
class(ucodones)
dim(ucodones)
rscuglobal=uco(unlist(s),index='rscu')
rscuglobal
sapply(s,uco(index='rscu')) #no se puede pasar parametros con definicienes dentro

rscu=function(seq){
  RSCU=uco(seq,index='rscu')
  return(RSCU)
} #se hace una funcion para meter en el sapply
usorelativo=sapply(s,rscu)
usorelativo
dim(usorelativo)
rownames(usorelativo)
