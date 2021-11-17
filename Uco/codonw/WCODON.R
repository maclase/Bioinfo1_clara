library(seqinr)
#los primeros ejes son los mar importantes, el coa reduce la tbala
#inicial de coordenadas, deja como primero al de mayor varianza
#que es el que tiene mas infoemacion
coa <- read.table('genes.coa',header = T)
genes <- read.table('CP003773.out', header = T)
codones <- read.table("codon.coa", header = T)


plot(codones$Axis1,codones$Axis2, pch='.')
text(codones$Axis1,codones$Axis2,label=codones$label, pos=4,cex=0.8)

bicho<- read.fasta('~/Bioinfo1_2021/Uco/codonw/CP003773.ffn')
largos <- sapply(bicho, length)

plot(coa$Axis2, coa$Axis2)
plot(coa$Axis1,coa$Axis4)
plot(coa$Axis1, largos)

plot(codones$Axis1, codones$Axis2)

cor(largos, x = coa$Axis1) #cuanto mas cercano a 1 mas lineal la correlacion
cor.test(largos, x = coa$Axis1)
cor(genes$Aromo, x = coa$Axis4)
plot(codones$Axis1, largos)

cor(largos, x =coa$Axis1)
