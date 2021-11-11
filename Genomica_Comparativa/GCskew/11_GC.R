

# Ejemplo vibrio de grafico oriloc ----------------------------------------


gbk = "vibrio.gbk"
bacdf = oriloc(gbk = gbk)
#bacdf <- oriloc(gbk="CP007025.gbk") # el nombre bacdf es por "bacteria data.frame". 
# GCA_000005845.2_ASM584v2_genomic.gbff
draw.oriloc(bacdf, main = "oriloc de vibrio")

#acumulacion se hace con cumsu
#Hagamos 
?draw.oriloc
bacdf
# todo eso se puede modificar:
draw.oriloc(ori, main = "Oriloc vibrio",
            xlab = "Map position in Kb",
            ylab = "Cumulated combined skew in Kb", las = 1, las.right = 3,
            ta.mtext = "Cumul. T-A skew", ta.col = "pink", ta.lwd = 1,
            cg.mtext = "Cumul. C-G skew", cg.col = "lightblue", cg.lwd = 1,
            cds.mtext = "Cumul. CDS skew", cds.col = "lightgreen", cds.lwd = 1,
            sk.col = "black", sk.lwd = 2,
            add.grid = TRUE)
draw.oriloc(bacdf,cg.col="blue",cg.lwd=2,sk.col="gray",sk.lwd=1) # Este es el comando que usaremos