p = rnorm(1000000,100,20)
hist(p,100)
mean(p)
sd(p)
# p = sample(1:20,1000000,replace=TRUE)
# p = sample(1:20,1000000,replace=TRUE,prob=c(21:40))

# De esa población p voy a hacer mil muestras de n=100
lista=list()
for(i in 1:1000){
  lista[[i]]=sample(p,100,replace=FALSE)#replace false porque son muestras reales, sin reposicion
}

medias = sapply(lista,mean)
# Ahora veo qué tanto varían las medias de las distintas muestras
sd(medias)
hist(medias,100)
# Pero si tengo una muestra sola, por ejemplo, la primera que hice

s1 = lista[[1]]

# En el Bootstrap lo que hago es un falso remuestreo: 
# tengo una única muestra original, de tamaño n, y hago muchas 
# submuestras con reposición, también de tamaño n. 
# Después, como hice antes (en las muchas muestras verdaderas), 
# le calculo la media a cada una, y después el sd de las medias.
# Acá calcularé las medias a medida que las voy haciendo
# e iré guardando solo eso, porque las muestras en sí no me interesan.
boot = function(vector,n){
  m=c()
  for(i in 1:n){
    v=sample(vector,length(vector),replace=TRUE) #en boot se hace con reposicion, muestras no reales
    m[i]=mean(v)
  }
  return(m)
}
mediasboot = boot(s1,1000)
sd(medias)
sd(mediasboot)
hist(mediasboot,100, col = 5)
par()
hist(medias,100, col = 3)

