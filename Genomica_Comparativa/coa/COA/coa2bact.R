library(seqinr)
library(ade4)

nombrearchivo1 <- "Listeria.cds.nc"
nombrearchivo2 <- "Pseudomonas.cds.nc"
bicho1<-read.fasta(nombrearchivo1)
bicho2<-read.fasta(nombrearchivo2)
table1GC<-data.frame(GC=sapply(bicho1,GC),GC1=sapply(bicho1,GC1),GC2=sapply(bicho1,GC2),GC3=sapply(bicho1,GC3))
table2GC<-data.frame(GC=sapply(bicho2,GC),GC1=sapply(bicho2,GC1),GC2=sapply(bicho2,GC2),GC3=sapply(bicho2,GC3))

gece1 <- GC(unlist(bicho1))
gece2 <- GC(unlist(bicho2))

dim(table1GC)
dim(table2GC)
tabGC <- rbind(table1GC,table2GC) #rbind= toma un vector, matriz o df y combina filas en este caso (r)
dim(tabGC)


mkdata <- function(seqs){
  tab <- sapply(seqs, uco)  
  tab <- as.data.frame(tab)
  return( tab )
}

tab1 <- mkdata(bicho1)
tab2 <- mkdata(bicho2)
dim(tab1)
dim(tab2)
tab <- cbind(tab1,tab2) #se pegan columnas aca
dim(tab) #64 codones=filas, 8439 genes=columnas

vector_b1 <- 1:ncol(tab1) # si bicho1 tiene 5000 genes, esto va de 1 a 5000.
vector_b2<- (ncol(tab1)+1):ncol(tab) # si bicho 2 tiene 2000 genes, esto va de 5001 a 7000.

coa_ <- dudi.coa(tab, scan = FALSE, nf = 3)
facaa <- as.factor(aaa(translate(sapply(rownames(tab),s2c))))
scua_ <- wca(coa_, facaa, scan = FALSE, nf = 3)
scua <- scua_
attributes(scua_)

# Esto que sigue es para automatizar el xlim y el ylim, y no tener
# que ponerlos a mano cada vez que cambiamos de bacterias.

xmax_ <- max(scua_$co[,1])
xmin_ <- min(scua_$co[,1])
ymax_ <- max(scua_$co[,2])
ymin_ <- min(scua_$co[,2])
xmin <- xmin_
ymin <- ymin_
xmax <- xmax_
ymax <- ymax_
# scua$co es una dataframe con los nombres de los genes (de ambas especies), 
# y las coordenadas de los puntos en los tres primeros ejes.
# scua$li es lo mismo pero con los nombres y coordenadas de los codones en los mismos tres ejes
dev.new()
pdf("coa_2_bacterias.pdf")
plot(scua_$co[vector_b1,1],scua_$co[vector_b1,2],pch=20,col="orange",xlab="Axis 1",ylab="Axis 2",
    xlim=c(xmin,xmax),ylim=c(ymin,ymax))
par(new=TRUE)
plot(scua_$co[vector_b2,1],scua_$co[vector_b2,2],pch=20,col="sky blue",xlab="Axis 1",ylab="Axis 2",
    xlim=c(xmin,xmax),ylim=c(ymin,ymax))
s.label(scua_$li, add.plot = TRUE, clab = 0.75,boxes = FALSE,xlim=c(xmin,xmax),ylim=c(ymin,ymax))
abline(v=0,h=0,lty="dotted")
dev.off()

cor(tabGC,scua$co)

plot(coa_$co[vector_b1,1],coa_$co[vector_b1,2],pch=20,col="orange",xlab="Axis 1",ylab="Axis 2",
     xlim=c(xmin,xmax),ylim=c(ymin,ymax)) #agarra filas del bicho 1 =vector_b1
par(new=TRUE) #para generar dos graficos en una misma imagen, pero hay q definir limites de ejes, para que queden visibles los dos
plot(coa_$co[vector_b2,1],coa_$co[vector_b2,2],pch=20,col="sky blue",xlab="Axis 1",ylab="Axis 2",
     xlim=c(xmin,xmax),ylim=c(ymin,ymax)) #
s.label(coa_$li,add.plot=TRUE,clab=0.75,boxes=FALSE,xlim=c(xmin,xmax),ylim=c(ymin,ymax))

