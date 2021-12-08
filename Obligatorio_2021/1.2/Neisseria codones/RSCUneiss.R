
library(seqinr)

fasta <- read.fasta("NeisseriaCDS.ffn") #Leemos el archivo fasta

fasta_modulo3 <- fasta[-which(getLength(fasta)%%3 >0)]
length(fasta) # El objeto fasta tiene 4045 secuencias (CDSs)

sec <- unlist(fasta_modulo3) # crea un elemento unico, desglozando la lista original
class(sec)
length(sec)

#uco cuenta las veces que aparecen los codones (nucleotidos de a tres)  en una secuencia
suco <- uco(sec)

#Se procesa lo obtenido con uco
#creeando un data frame en donde se realacionen las frecuencias de aparicion de los 
#codones y los aminoacidos con sus nombres
freqaa <- as.data.frame(suco)
#Se agregan nombres de aminoacidos
nombresaa  <-  translate(sapply(names(suco),s2c))
freqaa[,3] <- nombresaa
freqaa[,4] <- aaa(freqaa[,3]) 
#Se ordenan las columnas teniendo en cuenta la frecuencia de aparicion
freqaa_ord <- freqaa[order(freqaa[,2]),]
freqaa_ord
#Se determina el aa de mayor frecuencia
codon_mas_freq <- freqaa[which.max(freqaa$Freq), 1]
codon_mas_freq
#En este caso el codon de mayor frecuencia es el gaa (Glu)


#Se grafican las frecuencias de los codones en un grafico de barras
barplot(freqaa$Freq, names.arg = df$V3, cex.names = 0.6, 
        main = "Frecuencia de uso de codones Nesseria", col = 4, ylab = "Frecuencia", xlab = "Codones")
