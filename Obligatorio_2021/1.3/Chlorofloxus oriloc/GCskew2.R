library(seqinr)

gbk = "Chloroflexus.gbk"
chlordf= oriloc(gbk = gbk)
draw.oriloc(chlordf)

?draw.oriloc
chlordf

draw.oriloc(ori = chlordf, main = "Sesgo Chloroflexus", xlab = "posicion genomica (Kb)",
            ylab = "sesgo acumulado (Kb)",las = 1, las.right = 1, ta.mtext = "TA acumulado",
            ta.col = "red", cg.mtext = "GC acumulado", cg.col = "blue", cds.mtext = "CDS acumulado",
            cds.col = "green", sk.col = "black", sk.lwd = 2,add.grid = T)
#El origen de replicacion se ve aproximadamente en 2000 Kb, donde coincen el punto bajo de el sesgo TA y CDS,
#junto con el maximo de sesgo GC.La terminacion de la replicacion se puede estimar en las 4500 Kb.