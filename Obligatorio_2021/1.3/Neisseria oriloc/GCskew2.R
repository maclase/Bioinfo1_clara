library(seqinr)

gbk = "Neisseria.gbk"
neissdf= oriloc(gbk = gbk)
draw.oriloc(neissdf)

?draw.oriloc
neissdf

draw.oriloc(ori = neissdf, main = "Sesgo Neisseria", xlab = "posicion genomica (Kb)",
            ylab = "sesgo acumulado (Kb)",las = 1, las.right = 1, ta.mtext = "TA acumulado",
            ta.col = "red", cg.mtext = "GC acumulado", cg.col = "blue", cds.mtext = "CDS acumulado",
            cds.col = "green", sk.col = "black", sk.lwd = 2,add.grid = T)

#La terminacion de replicacion se ve aproximadamente en 750 Kb, donde coincen el punto alto de el sesgo TA y CDS,
#junto con el minimo de sesgo GC. El inicio de la replicacion se puede estimar en las 1800 Kb.