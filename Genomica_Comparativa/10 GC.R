
# Crear una secuencia aleatoria para el calculo de sesgo C/G --------------
library(seqinr)
seq <- sample(x=c('a','c','g','t'),size = 100, replace = T) #replace=TRUE para que sea con reposicion

# Calcular composiciones nucleotidicas de ventanas ------------------------

f1 <- seq[1:10]
f2 <- seq[11:20]
f3 <- seq[21:30]
f4 <- seq[31:40]
f5 <- seq[41:50]

lista <- list(f1,f2,f3,f4,f5)
lista
contar_bases <- function(secuencia) {
  c <- count(seq=secuencia, wordsize = 1,alphabet = s2c('acgt')) #hasta aca te dice cuanto hay de cada base
  csobreg <- conteos[2]/conteos[3] #aca hace la division C/G
  return(csobreg)
}
conteos <- count(seq=lista[[1]], wordsize = 1,alphabet = s2c('actg'))
contar_bases(f1) 

sapply(lista, contar_bases)

# Con una secuencia mas larga ---------------------------------------------
#A cada base le asiganamos un numero del 1 al 100, para definir a que ventana pertenece cada base
seq2 <- sample(x=c('a','c','g','t'),size = 10000, replace = T) 
count(seq2, wordsize = 1)
wsize <- 100 #tamano de las ventanas
len <- length(seq2)/wsize

#Hacemos 100 ventanas de 100 bases
#a cada base le asignamos un numero del 1 al 100, para definir a que ventana pertenece

f<- list() #lista vacia
for (i in 1: round (length(seq2)/wsize)){ #round por si queda decimal el resultado de la division
  f[[i]]=rep(i,wsize)
} 
f
f <- unlist(f) #unlist desarma la lista y se concatenan todos los elementos en un vector

#split separa cada elemento en una lista, segun el argumento en f
ventanas <- split(x=seq2, f= f) #obtendo al secuencia de cada ventana
contar_bases <- function(secuencia){
  conteos <- count(seq=secuencia, wordsize = 1,alphabet = s2c('actg')) #hasta aca te dice cuanto hay de cada base
  #csobreg <- conteos[2]/conteos[3] #aca hace la division C/G
  cmenosg <-conteos[2]-conteos[3]
  gc_skew <- (conteos[3]-conteos[2])/(conteos[3]+conteos[2])
  return(gc_skew)
}
CsobreG <- sapply(ventanas, contar_bases) #contar bases es la funcion de arriba, debuelve la divicion C/G para cada ventana
CsobreG
CmenosG<- sapply(ventanas, contar_bases)
GCskew <- sapply(ventanas, contar_bases)
GCskew
# Plotear C/G -------------------------------------------------------------

plot(CsobreG)
plot(CmenosG, type = 'l')
plot(GCskew)
# Determinacion del sesgo C/G de un genoma --------------------------------

genoma <- read.fasta('Ypestis.fna')

#definimos largo de las ventanas
getLength(genoma)
wsize <- 10000

#se crean ventanas 

f <- list()
for (i in 1:round(length(genoma)/wsize) ) {
  f [[i]]=rep(i,wsize)
}
f
f <- unlist(f)
ventanas <- split(x=seq2, f= f) #obtendo al secuencia de cada ventana

CsobreG <- sapply(ventanas, contar_bases) #contar bases es la funcion de arriba, debuelve la divicion C/G para cada ventana
CsobreG


plot(CsobreG, type='1')
