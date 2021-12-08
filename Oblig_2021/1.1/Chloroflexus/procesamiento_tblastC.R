library(seqinr)
tablast_Chloro <- read.table("tablast_Chloro_argC", sep="\t") #Se lee el archivo generado por BLASTP previamente 
#Se genera un vector para nombrar las columnas
colnames(tablast_Chloro) <- c("query","subject","porcid","len_aln","mism","gapop",
                             "qstart","qend","sstart","send","evalue","score", "qlen", "slen"
)
head(tablast_Chloro)

#Se filtra la tabla segun dos criterios seleccionados para determinar presencia de 
#secuencia de resistencia a antibioticos: 
# 1) El porcentaje identidad>50% 
tC_50 <-tablast_Chloro[which(tablast_Chloro$porcid>50), ] 

#2) Largo de alimiento por lo menos la mitad del largo del subject
tC_len <- tablast_Chloro[which((tablast_Chloro$slen*0.5)<(tablast_Chloro$len_aln)), ]

#Se plantea filtro segun las dos condiciones
tC_50_len <- tC_50[which((tC_50$slen*0.5)<(tC_50$len_aln)), ] 
#en vez de desarrollar lo del parentesis se podria usar tC_len

#guarda el resultado del filtrado en un archivo
write.table(x = tC_50_len, file = "tablast_Chloro_argC_filtrado.tsv", quote = F, 
            sep = '\t', row.names = F)
#Se opta por guardar la tabla filtrada por las dos condiciones auqnue tambien se podria
#generar archivos de estas por separado con las variables tC_50 y tC_len


# RESULTADOS --------------------------------------------------------------

#Se encontro 1 cds que coincide con las secuencias reportadas que infieren resistencia a 
#antibioticos a las bacterias segun el archivo arg.fasta
#El porcentaje de identidad es de 51,8%