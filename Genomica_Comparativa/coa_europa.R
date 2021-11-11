 # CASO DE ANÁLISIS DE CORRESPONDENCIA (COA) NO BIOLÓGICO, PARA ENTENDER LA IDEA. 
# (Ejemplo modificado, tomado de internet, ni sé si
# es una encuesta real o un ejemplo inventado para enseñar estadística)

# Tabla: opinión de ingleses sobre habitantes de distintos países (incuyendo Inglaterra).
# No sé cómo era, pero supongamos que de una lista de características personales, 
# les hacían elegir las tres que considerara más descriptivas de la gente de cada país.
# En la tabla original, las filas corresponden a países, y las columnas a características. 
# Los valores indican cuántas veces se mencionó esa característica asociada con ese país.
# Lo imnportante es que de la tabla original uno puede, con suerte, deducir algunas cosas, 
# pero después de hacer el COA todo queda más claro.

library(ade4)
t<-read.table("tabla_europa",header=T)
# Observamos la tabla
t

# Hacemos el análisis de correspondencia
coa <- dudi.coa(t)
2

# Dejo que me muestre el barplot de valores propios, que indica visualmente la importancia de 
# cada eje como explicador del sistema. Cuando me pregunta, y tras mirar el barplot, le digo 
# que me de solo dos ejes. Veremos que hay otra forma de hacer esto.
# EN LOS ESTUDIOS QUE HACEMOS NOSOTROS, CADA PAÍS (COLUMNA) SERÍA UN GEN, 
# Y CADA CARACTERÍSTICA (FILA), UN CODÓN, Y LOS VALORES SERÍAN,
# POR EJEMPLO, LA FRECUENCIA DEL CODÓN DE ESA FILA EN EL GEN DE ESA COLUMNA.

attributes(coa)

coa$li  # coordenadas de cada país en los ejes 1 y 2, en el mismo gráfico. GENES, PROTEÍNAS, ETC.
coa$co  # coordenadas de cada "característica" en los ejes 1 y 2. CODONES, AMINOÁCIDOS, ETC.

# Haré un plot de las "características"
eje1 <- coa$co$Comp1
eje2 <- coa$co$Comp2

plot(eje1,eje2,pch=20,xlim=c(-1,1),ylim=c(-1,1),xlab="Eje 1",ylab="Eje 2") 
# Cada punto corresponde a una característica (valiente, eficaz, vago). En lo nuestro: cada punto, un gen.

abline(h=0,v=0,col="lightblue")

# El gráfico, como está, no nos dice mucho, pero agreguemos los labels por PAÍS (serían las etiquetas
# de los codones: AAA,AAC,AAG...)

s.label(coa$li, add.plot = TRUE, clab = 0.75,boxes = TRUE)

# Acá ya podemos encontrar algo importante: en el eje1, del lado derecho, están los países latinos, 
# y del izquierdo, los no latinos.
# En el eje2, se ubican arriba los países exitosos, o económicamente poderosos. 
# Cabe destacar que los datos son de una época en que España estaba notoriamente peor que Italia 
# (los jugadores de fútbol más cotizados, por ejemplo, jugaban el la liga italiana).

# Aquí, en que la tabla es pequeña, podemos agregar, además, las etiquetas de las características. 
s.label(coa$co, add.plot = TRUE, clab = 0.75,boxes = FALSE)

# En bioinfo, eso llenaría  el gráfico de nombres de genes (lo cual sería un lío), 
# pero si tuviéramos interés en estudiar determinado gen o grupo de ellos, podríamos
# ubicarlo solo a él en el gráfico (me refiero a las etiquetas), conociendo su ubicación en el fasta original.
# Por ej, si son los tres primeros genes del fasta, pondríamos s.label(coa$li[1:3], etc.).
# También podemos "pintar" los puntos.
# Supongamos que lo que acá son las características más asociales a las estrellas de cine (elegante, arrogante, sexy) correspondiera,
# en un estudio bioinformático con miles de puntos, a los "genes de alta expresión". Los destacamos con color.


plot(eje1[1:3],eje2[1:3],pch=20,col="red",xlim=c(-1,1),ylim=c(-1,1),xlab="Eje 1",ylab="Eje 2")
abline(h=0,v=0,col="lightblue")
par(new=TRUE)
plot(eje1[-c(1:3)],eje2[-c(1:3)],pch=20,xlim=c(-1,1),ylim=c(-1,1),axes=FALSE,ann=FALSE) 
# no redibujo ni reanoto los ejes

# Volvemos a etiquetar
s.label(coa$li, add.plot = TRUE, clab = 0.75,boxes = TRUE)
s.label(coa$co, add.plot = TRUE, clab = 0.75,boxes = FALSE)

# Es notable que Inglaterra está un poco más abajo (en el eje2) de lo que uno la ubicaría; 
# probablemente sea un tema de autopercepción; a todos nos gusta quejarnos de lo mal que estamos.
# La distancia de cada "característica" con el país correspondiente indica cuánto hay de esa 
# característica en ese país, siempre según el imaginario colectivo de los encuestados.

# En nuestros análisis, los ejes se suelen vincular a características como el contenido GC, 
# la hidrofobicidad de los aminoácidos, etc.

# hay dos opciones: scan = FALSE, nf = 3 (u otro número). La primera evita el plot inicial de los eigenvalues, 
# y la segunda dice, de entrada, cuántos ejes querés para hacer los plots. 
# Lo usual, cuando se hace el coa o el wca con uso de codones, son tres ejes.


