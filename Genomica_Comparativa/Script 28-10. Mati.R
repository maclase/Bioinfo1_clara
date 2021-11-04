library(seqinr)
seq = s2c("actacgatcgatagctagctacgataatgcgcgctagctagatgatattcagagtctactaggggact") # Tenemos que pasar a con la funcion s2c a strings
seq2 = sample(s2c("actg"), size = 10000, replace = TRUE)


Fun = function(secuencia){
  contar = count(secuencia,wordsize = 1,by = 1)
  fraccion = contar[2] / contar[3]
return(fraccion)
  }
Fun (seq2)


# Utilizamos el script de composicional2 ----------------------------------

length(seq2)

tamaño_ventanas=100 # Para estudiar el contenido gc, lo hago por ventanas.
cantidad_ventanas = length(seq2)/tamaño_ventanas # 100 aporx de ventanas esta bien
# Ventanas son fragmentos del genoma para estudiarlo de a poco.


f=list()
for (i in 1: round(length(seq2)/tamaño_ventanas)) { 
  f[[i]] = rep(i,tamaño_ventanas)
} 
# El round esta por si me queda un numero decimal, para redondear al entero mas cercano si me queda una ventana con coma
# Voy diciendole a cada base a que numero de ventana pertenece. Creo una lista donde me dice cada base a que ventana pertenece.
#Loop para numeros grandes, divido en segmentos.

length(f)
f=unlist(f) # Creamos un unico string, desligamos la lista
length(unique(f))
length(f)

ventanas = split(x = seq2, f = f) # Con la funcion split lo que hacemos es por un lado definimos la secuencia que estamos trabajando (seq2), y a partir de esto distribuimos las bases de la secuencia en base a
# En base al vector que habiamos contruido con el loop y unlist donde tendremos 111111111111222222222222222222223333333333333334444444444444
# Me quedaran una lista con todas las ventanas con cada una sus bases.

# En el script habiamos hecho esto, ahora cambiamos 
#sp=list()
#sp=split(g,f) # Función SLPIT sirve para dividir todo.
#length(sp) # Tenemos 93 ventanas

C_sobre_G = sapply(ventanas, Fun) # Por medio de la funcion sapply aplico la funcion que construimos (F) a la lista de las ventanas que fabrique
plot(C_sobre_G)
plot(C_sobre_G, type = "h")
barplot(C_sobre_G)
