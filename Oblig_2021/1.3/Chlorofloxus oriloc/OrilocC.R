library(seqinr)

#Previamente se obtiene el archivo gbk descomprimiendo el archivo descargado
#Se abre ek archivo de gene bank y se ejecuta oriloc
archivo_gbk = "Chloroflexus.gbk"
chlordf= oriloc(gbk = archivo_gbk)

chlordf
#Se grafica los sesgos de GC y AT
draw.oriloc(ori = chlordf, main = "Sesgo Chloroflexus", xlab = "posicion genomica (Kb)",
            ylab = "sesgo acumulado (Kb)",las = 1, las.right = 1, ta.mtext = "TA acumulado",
            ta.col = "red", cg.mtext = "GC acumulado", cg.col = "blue", cds.mtext = "CDS acumulado",
            cds.col = "green", sk.col = "black", sk.lwd = 2,add.grid = T)
#El origen de replicacion se ve aproximadamente en 2000 Kb, donde coincen el punto bajo de el sesgo TA y CDS,
#junto con el maximo de sesgo GC.La terminacion de la replicacion se puede estimar en las 4500 Kb.