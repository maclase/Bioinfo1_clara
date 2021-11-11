library(seqinr)

# El archivo bicho1.ffn contiende CDS de una bacteria Veamos si  --------
#  los codones sinónimos son usados con frecuencias similares. -----------

fasta=read.fasta("bicho1.ffn")
length(fasta)

f10=fasta[[10]]
f10
length(f10)/3 #para saber si es sec codificante, da la cantidad supuesta de aa
cuenta10=count(f10,wordsize = 3,by = 3) #para contar codones, count de 3 en 3, no se especifica alfabeto
cuenta10


# uco hace lo mismo que count con wordsize=3 y by=3 ----------------------

uco10=uco(f10)
cuenta10;uco10
identical(cuenta10,uco10) #pregunta si el lo mismo lo obtenido con count y con uco


# Voy a ver a qué aminoácido corresponde cada codón -----------------------
names(uco10)
namesaa = translate(sapply(names(uco10),s2c)) #s2c transforma en vectores cortos
dfv=as.data.frame(uco10) #hacer data frame
dfv[,3]=namesaa #se agrega como columna 3 los nombres de aa q creamos
dfv
dfv=dfv[order(dfv[,3]),] #Se ordena tabla por la columna 3 
dfv
barplot(dfv$Freq,names.arg = dfv$V3,cex.names = 0.6, col = 4, main = "bicho")


# hagamos lo mismo para toda la secuencia (UCS global)

sec=unlist(fasta) #se unen todos los vectores (secuencias) de una lista
ucosec=uco(sec)
df=as.data.frame(ucosec)
df[,3]=namesaa
df
df=df[order(df[,3]),]
barplot(df$Freq,names.arg = df$V3,cex.names = 0.6, col = 3, main = "bicho")


# Relative Synonimous Codon Usage RSCU ------------------------------------

# En este caso se calculan las frecuencias pero dentro de los sinónimos
# de cada aminoácido. Después se normaliza para que la suma total de los
# RSCU de un aá sea igual a la cantidad de sinónimos del aá. 

rscu=uco(sec,index = "rscu")
rscu
sum(rscu)
df=as.data.frame(rscu)
df[,2]=namesaa
df
df=df[order(df[,2]),]
barplot(df$rscu,names.arg = df$V2,cex.names=0.6, main = "rscu de bicho", col = 7) #columna mas alta sigue siendo L, 
                                                 #la cantidad q hay es por un sesgo importante, no tanto por ser un aa comun
df$V2
col=c(1,1,1,1,2,2,2,2,2,2,3,3,4,4,5,
      5,6,6,7,7,8,8,1,1,1,1,2,2,3,3,3,4,4,4,4,4,4,5,5,6,7,7,8,8,8,8,
      1,1,1,1,1,1,2,2,2,3,3,3,3,4,5,5,6,6,6,6)
length(col) 
