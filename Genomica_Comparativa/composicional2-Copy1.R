library(seqinr)
g=read.fasta("Ypestis.fna")[[1]] #lo convertimos en vactor para poder aplicar algunas funciones
gc=GC(g)
gc


# Ejemplo de hacer el dataframe para dividir secuencia y poder calcular para cada segmento ------------------------------------------

wlen=50000 #largo de ventanas, si es 500000 es muy bajo, si pongo 500 es muy alto
length(g)
length(g)/wlen
letras=letters
length(letras)
vector=s2c('11122233344455566677788899')
length(vector)
split(letras,vector)
data.frame(letras,vector)

#para fabricar numero de ventanas q permitan dividir la secuencia automaticamente, sin ingresar el vector
f=list()
for (i in 1: round(length(g)/wlen) ){
  f[[i]] = rep(i,wlen)
}
length(f)
f=unlist(f)

length(unique(f))
length(f)
sp=list()
sp=split(g,f) #sp es el cromosoma dividido en partes iguales
length(sp)

gcsp=sapply(sp,GC)
gcsp
summary(gcsp)
plot(gcsp,type="l") #no me corre dice q los margenes son muy grandes



#otro modo
plot(gcsp,pch=20,col="red")
par(new=TRUE)
plot(gcsp,type="h")

#hagamos un vector de colores
colores=rep("black",length(sp))
alto=grep(TRUE,gcsp>0.49) #le da a true a los q cumplen la condicion dada
colores[alto]="red"
medio=grep(TRUE,gcsp<=0.49 & gcsp > 0.46)
colores[medio] = "orange"
bajo=grep(TRUE,gcsp<=0.46)
colores[bajo]="light blue"
#al usar el vector colores como parametro para los colores
#del plot, quedan coloreados diferencialmente
plot(gcsp,pch=20,col=colores)
la par(new=TRUE)
plot(gcsp,type="h")



# hasta acá








cmenosg=function(seqdna){
  cuenta= count(seqdna,wordsize = 1,start = 0,by = 1,alphabet = s2c("acgt"))
  return(cuenta[2]-cuenta[3])
}

cmg=sapply(sp,cmenosg)
mean(cmg)
plot(cmg,type="l")
abline(h=mean(cmg))


cds=read.fasta("Ypestis.cds")
cds=unlist(g)
length(cds)
gccds=GC(cds)
gccds
fcds=list()
for (i in 1: round(length(cds)/wlen) ){
  fcds[[i]] = rep(i,wlen)
}
spcds=split(cds,fcds)
cmgcds=sapply(spcds,cmenosg)
mean(cmgcds)
mean(cmg)
