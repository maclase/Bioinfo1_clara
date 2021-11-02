

# #Ejemplo hecho por nosotros ------------------------------------ --------



largos_ej<- c(5.12, 3100, 0.81, 3.12, 6.44, 1.7, 1.7, 472, 143.57, 4.53, 4.77, 32000, 2226, 2037)
names(largos_ej) <- c("E.coli", "Sapiens", "M. pneumo", "rodo", "babesia", "marinus","mathanococcus","phys","droso","methanosarcina","helecho", "salamandra", "garrapata", "blatella")
barplot(largos_ej,names.arg = names(largos),space = 2)


# Aca empieza ejemplo de clase, con tabla de Bioinfo1_2021 ----------------



largos <- read.table("LargosGenomas.txt", sep="\t", stringsAsFactors = F)
barplot(largos$V2, names.arg=largos$V1, log="y", col = palette())

##Buscar bacterias unicamente, con which buscamos las q no tienen nad aen columna 3

which(largos$V3=="")
largos_bacterias <- largos[which(largos$V3!=""),] ##aca quedaron guardados los genomas de bacterias unicamente, filtrando segun tengan o no info en tercera columna
largos_bacterias<- largos_bacterias[-8,]
dim(largos_bacterias)
barplot(largos_bacterias$V2, names.arg = largos_bacterias$V1, col = c(1,2,3,4,5), cex.names = 0.6)




