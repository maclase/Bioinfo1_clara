library(seqinr)

cr=read.fasta("ChrNitrobacter.ffn") #cr=lista
p1=read.fasta("Pls1Nitrobacter.ffn") #p por plasmido
p2=read.fasta("Pls2Nitrobacter.ffn")
p3=read.fasta("Pls3Nitrobacter.ffn")

length(cr)
length(p1)
length(p2)
length(p3)

crlen=sapply(cr,length)
mean(crlen)
p1len=sapply(p1,length)
mean(p1len)
p2len=sapply(p2,length)
mean(p2len)
p3len=sapply(p3,length)
mean(p3len)
listalen=list("crlen"=crlen,"p1len"=p1len,"p2len"=p2len,"p3len"=p3len) #lista de vectores de diferente largo, no puede ser matriz  dataframe
boxplot(listalen,outline = FALSE,col=c(2,3,4,5)) #si se aplica boxplot a los elementos de la lista grafica las cajas de c/u en una sola grafica

crgc1 = sapply(cr,GC1)
crgc2 = sapply(cr,GC2)
crgc3 = sapply(cr,GC3)
p1gc1 = sapply(p1,GC1)
p1gc2 = sapply(p1,GC2)
p1gc3 = sapply(p1,GC3)
p2gc1 = sapply(p2,GC1)
p2gc2 = sapply(p2,GC2)
p2gc3 = sapply(p2,GC3)
p3gc1 = sapply(p3,GC1)
p3gc2 = sapply(p3,GC2)
p3gc3 = sapply(p3,GC3)
listagc = list("crgc1"=crgc1,"crgc2"=crgc2,"crgc3"=crgc3,
               "p1gc1"=p1gc1,"p1gc2"=p1gc2,"p1gc3"=p1gc3,
               "p2gc1"=p2gc1,"p2gc2"=p2gc2,"p2gc3"=p2gc3,
               "p3gc1"=p3gc1,"p3gc2"=p3gc2,"p3gc3"=p3gc3) #al igual q arriba se hace una lista con todos los GC
boxplot(listagc,col=c(2,2,2,3,3,3,4,4,4,5,5,5)) #todos los GC (1,2,3) de un mismo replicon con el mismo color, por eso repite
boxplot(listagc,col=c(2,2,2,3,3,3,4,4,4,5,5,5),outline = F) #le saco los GC mas lejanos de la medio (caja), porque no aportan

# aminoácidos

crpep=sapply(cr,translate)
p1pep=sapply(p1,translate)
p2pep=sapply(p2,translate)
p3pep=sapply(p3,translate)
secpep=crpep 
contaraa=function(secpep){
  cuenta=count(seq = secpep,wordsize = 1,start = 0,by = 1,freq = TRUE,alphabet = a()[-1]) #
  return(cuenta)
}
MCcr=t(sapply(crpep,contaraa))
class(MCcr) #MATRIZ
round(head(MCcr),2) # esto es solo para ver más claro, pero dejo en la matriz los valores más exactos
MCp1=t(sapply(p1pep,contaraa))
MCp2=t(sapply(p2pep,contaraa))
MCp3=t(sapply(p3pep,contaraa))

boxplot(MCcr,main="MCcr",ylim=c(0,0.35))

boxplot(MCcr,outline=FALSE, main="MCcr",ylim=c(0,0.21))
boxplot(MCp1,outline=FALSE,main="MCp1",ylim=c(0,0.21))
boxplot(MCp2,outline=FALSE,main="MCp2",ylim=c(0,0.21))
boxplot(MCp3,outline=FALSE,main="MCp3",ylim=c(0,0.21))

p=as.vector(table(unlist(cr)))
lista=list()
for(i in 1:1000){
  lista[[i]]=sample(s2c("acgt"),size = mean(crlen)*sample(seq(0.6,1.4,by=0.1),1),
                    replace = TRUE,prob=p)
}

listapep=sapply(lista,translate)
Mlista=t(sapply(listapep,contaraa))
boxplot(Mlista,outline=FALSE,main="MCp3",ylim=c(0,0.21))
cor.test(colMeans(Mlista),colMeans(MCcr))
colnames(Mlista)
