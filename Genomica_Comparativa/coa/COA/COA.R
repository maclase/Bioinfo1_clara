library(seqinr)
library(ade4)

bicho<-read.fasta("bicho.ffn")

#se hace tabla GC, 1, 2 y 3. Tantas filas como genes y 4 columnas
tableGC<-data.frame(GC=sapply(bicho,GC),GC1=sapply(bicho,GC1),GC2=sapply(bicho,GC2),GC3=sapply(bicho,GC3)) 
head(tableGC)

# funcion para facilitar el uco de cada gen con sapply despues ------------

mkdata <- function(seqs){
  # seqs es una lista de seqs como la que se genera con read.fasta()
  tab <- sapply(seqs, uco)  
  tab <- as.data.frame(tab)
  return( tab )
}

tab <- mkdata(bicho) #generado con la funcion
dim(tab) #pàra hacer el coa codones=filas


# Se hace el COA y SCUA (coa con rscu) ------------------------------------

coa <- dudi.coa(tab, scan = FALSE, nf = 3) #dudi.coa es comando que usa el paquete ade4 para hacer coa
facaa <- as.factor((translate(sapply(rownames(tab),s2c)))) #traducir los nombres de tab a aminoacidos, primero pasar a vector, despues traducir
scua <- wca(coa, facaa, scan = FALSE, nf = 3) #wca = within coa -> hace lo que el rscu, sobre el coa directamente, sin chi cuadrado

# Datos que van a ir en los ejes del grafico despues, coordenadas en dataframes  --------

attributes(scua)

#scua$co -> coordenadas de los genes en el grafico, scua$co es data frame 
class(scua$co)
scua$co
dim(scua$co) # 2153 genes x 3 columnas

#scua$li -> coordenadas de los codones en el grafico, scua$li es data frame 
class(scua$li)
scua$li
dim(scua$li) #64 codones x 3 coordenadas


# realacionar los gc (1,2,3) con coordenadas de los genes en lso ejes -------------

round(cor(scua$co,tableGC),2) #round deja numeros redondeados segun las cifras que pases (2 en este caso)
plot(scua$co$Comp1, tableGC$GC3)


par(mfrow=c(1,2))
?par

# grafico todas las coordenadas de los genes en 3D------------------------------
#uso dos graficos porque tengo 3 coordenadas, grafico relacionandolas de a dos, coord 1 con 2 y depsues 1 con 3
plot(scua$co[,1],scua$co[,2],pch=20,col="orange",xlab="Axis 1",ylab="Axis 2")
s.label(scua$li, add.plot = TRUE, clab = 0.75,boxes = FALSE) #etiqueto los codones
abline(v=0,h=0)

plot(scua$co[,1],scua$co[,3],pch=20,col="orange",xlab="Axis 1",ylab="Axis 3")
s.label(scua$li, add.plot = TRUE, clab = 0.75,boxes = FALSE)
abline(v=0,h=0)

##################################
# Ahora vamos a ver qué pasa con alguna otra variable (para ver si los ejes 2 y 3 se correlacionan con algo)

# codonw bicho.ffn
# Después le digo que sí (ENTER) a lo que me pregunte.
# Elijo el menú 4, y dentro de él, las opciones gc, gc3 e hidropaticidad.
# Aprieto ENTER y luego R, para que corra el programa (hay que ir viendo qué opciones me da).
# Después ENTER otra vez, y q, para salir.
# Hago un head del archivo *out, y con read.table lo importo a R.
# Elimino la primera columna de out y guardo el resultado en out2.
# Recordemos que los tres ejes están en scua$co.
# Hagamos las correlaciones entre todas esas variables.
# Después creo otras variables, por ejemplo, largos. ¿Cómo se hace?
# Me fijo si la posición en el genoma tiene algo que ver con el GC, o con el largo.
# 
