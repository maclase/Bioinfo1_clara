library(seqinr)
tablast_Neiss <- read.table("tablast_Neiss_argN", sep="\t") #Se lee el archivo generado por BLASTP previamente 
#Se genera un vector para nombrar las columnas
colnames(tablast_Neiss) <- c("query","subject","porcid","len_aln","mism","gapop",
               "qstart","qend","sstart","send","evalue","score", "qlen", "slen"
)
head(tablast_Neiss)

#Se filtra la tabla segun dos criterios seleccionados para determinar presencia de 
#secuencia de resistencia a antibioticos: 
# 1) El porcentaje identidad>50% 
tN_50 <-tablast_Neiss[which(tablast_Neiss$porcid>50), ] 

#2) Largo de alimiento por lo menos la mitad del largo del subject
tN_len <- tablast_Neiss[which((tablast_Neiss$slen*0.5)<(tablast_Neiss$len_aln)), ]

#Se plantea filtro segun las dos condiciones
tN_50_len <- tN_50[which((tN_50$slen*0.5)<(tN_50$len_aln)), ] 
#en vez de desarrollar lo del parentesis se podria usar tN_len

#guarda el resultado del filtrado en un archivo
write.table(x = tN_50_len, file = "tablast_Neiss_argN_filtrado.tsv", quote = F, 
            sep = '\t', row.names = F)
#Se opta por guardar la tabla filtrada por las dos condiciones auqnue tambien se podria
#generar archivos de estas por separado con las variables tN_50 y tN_len

# RESULTADOS --------------------------------------------------------------

#Se encontraron 6 cds que coinciden con las secuencias reportadas que infieren resistencia a 
#antibioticos a las bacterias segun el archivo arg.fasta
#Los porcentajes de identidad de estas son mayores a 99%
