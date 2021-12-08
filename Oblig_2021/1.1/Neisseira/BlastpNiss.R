
# Procesamiento fasta antibioticos ----------------------------------------

filename="argN" 

library(seqinr)
fasta <- read.fasta(paste(filename,"pep",sep = ".")) #leer archivo fasta de secuencia aminoacidos 
FASTA <- sapply(fasta, toupper) #pasamos secuencia a mayuscula
len <-length(FASTA) 
head(names(FASTA))
nombres <- names(FASTA) #se guardan los headers de los cds para luego hacer una data frame

#funcion para generar nuevos nombres, identificandolos con a y el numero del cds correspondiente
nombres_nuevos <- c()
for(i in 1:len){
  nombres_nuevos[i] <-paste("a",i,sep=".")
}
head(nombres_nuevos)

#se genera el fasta nuevo con la secuencia en mayuscula con los headers simplificados 
write.fasta(FASTA,names =nombres_nuevos,file.out <- paste(filename,"pepmay",sep=".")) 

#se crea data frame guardando la relacione entre nombres viejos y nuevos, generando un archivo
guardar_nombres<-data.frame(nombres_nuevos,nombres)
head(dfnames)
write.table(guardar_nombres,quote=FALSE,row.names = FALSE,
            file=paste("nombres_argN"),sep="\t")
#De esta manera se apronta el archivo fasta para luego generar la base de datos, de manera que quede 
#leible cuando se presenten los datos en una table 

# Procesamiento fasta de Neisseria ----------------------------------------

filename="NeisseriaCDS" 

library(seqinr)
fasta <- read.fasta(paste(filename,"cds",sep = ".")) 
FASTA <- sapply(fasta, toupper)
len <-length(FASTA) 
head(names(FASTA))
nombres <- names(FASTA) #se guardan los headers de los cds para luego hacer una data frame

#funcion para generar nuevos nombres, identificandolos con a y el numero del cds correspondiente
nombres_nuevos <- c()
for(i in 1:len){
  nombres_nuevos[i]= paste("n",i,sep=".")
}
head(nombres_nuevos)


#Se genera el fasta con los headers nuevos y en mayuscula
write.fasta(FASTA,names=nombres_nuevos,file.out = paste(filename,"cdsmay",sep=".")) 
#Se traduce la secuencia nucleotidica
fastapep <- lapply(FASTA,translate)
#Se guarda traduccion en fasta
write.fasta(fastapep,names = nombres_nuevos,file.out = paste(filename,"pepmay",sep="."))

#Genera archivo con data frame de los nombre nuevos y viejos
nombres_Neiss<-data.frame(nombres_nuevos,nombres)
head(nombres_Neiss)
write.table(nombres_Neiss,quote=FALSE,row.names = FALSE,
            file=paste("nombres_Neiss",sep = "."),sep="\t")

#Este paso se realiza con el mismo objetivo que el anterior, pero no se destina para la 
#creacion de bases de datos, sino que se utiliza el archivo fasta con aminoacidos en el blast



# BLASTP ------------------------------------------------------------------
#en terminal

#makeblastdb -in argN.pepmay -dbtype prot -out argN_db
#blastp -query NeisseriaCDS.pepmay -db argN_db -out tablast_Neiss_argN -outfmt "6 std  qlen slen"