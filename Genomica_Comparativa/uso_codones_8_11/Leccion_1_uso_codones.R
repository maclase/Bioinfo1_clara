library(seqinr)
s <- read.fasta("NC_008095.ffn")
#chechear q la seq este bien
largos <- sapply(s,length)
sum(largos%%3) #es divisible entre 3
s1 = s[[1]]
u <-uco(s1)
length(u)
S <- unlist(s)
ucoS <- uco(S)
n <- names(ucoS)
n1 <- sapply(n, s2c)
class(n1)
typeof(n1)
dim(n1)

nAA <- apply(X = n1, MARGIN =2 , translate) #margin=2 para hacer translate en columna
nAAA <- aaa(nAA)
length(nAAA)
#para saber posiciones de aa
df <- data.frame(nAAA)
grep("Arg",nAAA)
grep("Stp",nAAA)
grep("Leu",nAAA)
grep("Phe",nAAA)

names(ucoS) <- nAAA
ucoS
uco_codonS <- ucoS
names(uco_codonS)<-n

vectArg <- grep("Arg",nAAA)
vectArg
cArg <- uco_codonS[vectArg]
sum(cArg)
#proporcion en que c/codon esta representado
(cArg/sum(cArg))*100
barplot(cArg)
#RSCU -> (ni*xij)/sum(xij)
#n=6
#x=2253
#sumatoria de todas las secuencias=227090
RCSU_aga=6*2253/227090
uco(seq = S, index = 'rscu')

#suma rscu para un genoma es siemrpe 64
suma <- sum(uco(seq = S, index = 'rscu'))


vectIle <- grep("Ile",nAAA)
vectIle
cIle <- ucoS[vectIle]
sum(cIle)

vectTrp <- grep("Trp",nAAA)
vectTrp
cTrp <- ucoS[vectTrp]
sum(cTrp)

#RSCU -> (ni*xij)/sum(xij)
#n=6
#x=2253
#sumatoria de todas las secuencias=227090