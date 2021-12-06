t=read.table("tablast_Chlor_argC", sep="\t")
colnames(t)= c("qseqid","sseqid","qlen","slen","length","pident",
               "evalue","bitscore")
head(t)

t_50 <-t[which(t$pident>50)] #which cuales son las poscion que cumplen con la condicion
#me da q no hay columnas, no hay ningun pident mayor a 50

t_50_len <- t_50[which(t_50$slen*0.5)<(t_50$length), ] #transformar los numeros y posciciones en la tabla
write.table(x = t_50_len, file = "tablast_Neiss_argN_filtrado.tsv", quote = F, sep = '\t', row.names = F)

#t_f <- t[which((t$slen*0.5)<(t$length)), ]
#write.table(x = t_f, file = "tablast_Chlor_argC_filtrado_sololen.tsv", quote = F, sep = '\t', row.names = F)