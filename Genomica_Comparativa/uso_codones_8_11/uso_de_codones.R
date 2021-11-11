# Hagan una carpeta llamada codones y pongan en ella algunos .ffn que tengan por ahí.
# Vamos a ver algunas cosas que se pueden hacer para ver cómo un organismo, en este caso una bacteria, usa diferencialmente (o no) sus codones sinónimos para codificar proteínas. 

library(seqinr)

################################### repaso de vectores #######################
# si a un vector, por ejemplo,
vector = 1:20
#le hago
vector > 10
# o
vector == 3
# me devuelve TRUE por cada vez que la condición se cumple. (Es como si hicier un "if").

#Así
aaa()
aaa()=="Arg"

# Esto nos permite una forma de usar el grep de R, que es parecido pero distinto al grep de la shell de linux
vectArg = grep(TRUE,aaa()=="Arg")
vectArg
# En este caso, sin embargo, es más fácil hacer simplemente.
vectArg = grep("Arg",aaa())
# Si en vez de == estuviéramos buscando valores 
# que cumplen determinada condición (>=, <, u otras más complejas), ahí hay que usar el grep con "TRUE".


# Por otro lado, si yo tengo el vector (que ya viene definido en R)
letters
#y hago
letters[c(1,3,4,26)]
# Me devuelve un vector con los elementos 1,3,4 y 26.

# Todo esto era necesario para entender algunos pasos que daremos ahora.
##################################################################################
s=read.fasta("NC_008095.ffn")

#¿Qué nos da la función "uco"?
# Veamos
s1 = s[[1]]
s1
uco(s1)

# uco, usado así, nos da el conteo de cuántos codones hay de cada tipo en una secuencia. 
# No discrimina entre codones sinónimos y no sinónimos.
# Ahora hagamos el uco para todos los CDSs juntos, concatenados en un "supergén".
S=unlist(s)
ucoS = uco(S)
ucoS

# también podría haber hecho el de cada CDS, haciendo
ucoCadaCDS = lapply(s,uco)
head(ucoCadaCDS)

# Pero volvamos al uco global, ucoS.
ucoS

# Lo interesante es comparar los números dentro del mismo aminoácido, que no están influidos por las distintas 
# frecuencias en que aparecen los aá. 
# Por ahora, hagamos

namesaa = as.factor(aaa(translate(sapply(names(ucoS),s2c))))
# Esto que parece tan complicado es simplemente para tener una variable que consista en la traducción 
# a aminoácidos (en código de tres letras, p. ej. Lys) de los codones, ordenados alfabéticamente, 
# que es como me aparecen en la salida del comando uco.

namesaa
# ¿En qué posición están los STOP?
grep("Stp",namesaa)
# ¿Y las argininas?
vectArg = grep("Arg",namesaa) # las guardo en un vector
vectArg 
# copio ucos en la variable ucoSaa y después le cambio los nombres-codones por nombres-aminoácidos.
ucoSaa = ucoS
names(ucoSaa) = namesaa


# corroboro
ucoSaa
# veamos las argininas
vectArg
ucoSaa[vectArg]
# pero lo voy a usar en el otro (en ucoS, en vez de ucoSaa), para, ahora que sé que son de arginina, ver qué codón es cada uno.
ucoArg = ucoS[vectArg]
ucoArg
sum(ucoArg)
ucoArgp100 = 100*ucoArg/sum(ucoArg)
ucoArgp100
# para que no tenga tantos decimales (en un porcentaje, más de dos números después de la coma suele ser una exageración), hacemos

ucoArgp100 = round( 100*ucoArg/sum(ucoArg),2 )
ucoArgp100
barplot(ucoArgp100)
legend("topleft",legend="uso de codones de Arg",bty="n")
# Veamos qué pasa con la composición de bases total.


############################################

# Otra forma, en vez de cambiar los nombres, era agregar los de los aminácidos como otra columna de una data.frame.
dfu = data.frame(ucoS,namesaa)
head(dfu)
# Ahora hagamos esto para todos los aminoácidos.
listau = list()
for(i in aaa()){
  listau[[i]] = dfu[ grep(i,dfu[,3]),2 ]
  names(listau[[i]]) = dfu[ grep(i,dfu[,3]),1 ]
}
listau
# Entonces, si queremos saber el uso de los codones de Val, hacemos
listau$Val

# Todo esto es una introducción. Ahora veremos cómo contar los codones de modo que sus frecuencias
# estén expresadas en función del aminoácido al que pertenecen (el conteo es el mismo).
# El método más conocido se llama RSCU ()
# En él, las frecuencias se muestran de modo que si un AA tiene cuatro codones, la suma de sus frecuencias dará 4.
?uco
rscuS = uco(seq=S,index="rscu")
rscuS
sum(rscuS)

dfr = data.frame(as.table(rscuS),namesaa) 
# esto porque uco nos da una table, pero con la opción rscu, nos da un vector
listar = list()

for(i in aaa()){
  listar[[i]] = dfr[ grep(i,dfr[,3]),2 ]
  names(listar[[i]]) = dfr[ grep(i,dfr[,3]),1 ]
}
listar


sum(listar$Met)
sum(listar$Phe)
sum(listar$Ile)
sum(listar$Val)
sum(listar$Ser)
sum(unlist(listar))

# De este modo, la suma total da 64, pero no importa si en la secuencia hay muchas o pocas argininas, metioninas
# o lo que sea; la frecuencia de cada uno de sus codones no se verá afectada por ello.