
library(seqinr)
p53 = read.alignment("p53.aln",format="clustal") #en vez de leer el fasta se lee directo el alineamiento
d1=dist.alignment(p53,matrix="similarity") #se hacen matrices con las distancias entre secuencias de aa de cada uno
d2=dist.alignment(p53,matrix="identity")
round(d1,2)
round(d2,2)
as.matrix(round(d2,2))
class(d1)

library(ape) #metodo para 
t1 = bionj(d1) # bionj CREA UN ARBOL A PARTIR DEL ALINEAMIENTO
t2 = bionj(d2)
plot(t1,type = "unr")
plot(t2,type = "unr")
#pdf("t1.pdf")
#plot(t1,type = "unr")
#dev.off()

# Si tuviéramos el alineamiento en fasta, igual podríamos leerlo así:
read.alignment("p53.mafft.fas",format="fasta")

# otra forma de leer un alineamiento y usar otras funciones 
library(phangorn)
ph <- read.phyDat("p53.mafft.fas", type = "AA", format="fasta")
ph
# Este objeto no se muestra, solo se ve alguna información sobre él.
# antes leía formato clustal, pero ahora, si bien sigue en el help, no.

class(ph)
dph=dist.ml(ph,model = "Blosum62")
class(dph)
class(d1)
# O sea: la matriz de distancias es el mismo tipo de 
# objeto en seqinr y en phangorn
tphUPGMA= upgma(dph)
class(tphUPGMA)
plot(tphUPGMA)

tphNJ=bionj(dph)
plot(tphNJ,use.edge.length = TRUE)

# BOOTSTRAP
# hago lo mismo para que me quede ordenado
ph <- read.phyDat("p53.mafft.fas", type = "AA", format="fasta")
dmph= dist.ml(ph,model="Blosum62")
tree=NJ(dmph)

NJtrees = bootstrap.phyDat(ph,
                           FUN =function(x)NJ(dist.ml(x)),bs = 1000)
plot(tphNJ)
treeNJ=plotBS(tphNJ, NJtrees, "phylogram",bs.col = "red")
treeNJ
attributes(treeNJ)
attributes(NJtrees)
# otro método
plot(tree)
fit= pml(tree,ph) # fit es otro árbol
plot(fit)
fit <- optim.pml(fit, rearrangement="NNI")
plot(fit)
plot(fit,"unr")
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
treeBS=plotBS(fit$tree,bs,bs.col = "red")
