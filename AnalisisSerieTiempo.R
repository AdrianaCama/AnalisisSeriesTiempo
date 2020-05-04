################################################################################################
################################################################################################
################################################################################################
# LAT4017 Series de Tiempo
# Universidad de las Américas Puebla
#
# Proyecto final:
#   Pronóstico de los precios de gas natural en 
#   Estados Unidos utilizando series de tiempo
#
# Integrantes: 
#   Adriana Camarillo Durán 15574
#   Ariel Arturo Ortega Alegría 155804
#
# Profesora:
#   Dra. Daniela Cortés Toto
#
# Fecha: 
#   04 de mayo de 2020
################################################################################################
################################################################################################
################################################################################################
################################################################################################

### Librerías ##################################################################################
library(forecast)
library(tseries)
library(fpp2)
library(ggplot2)
library(fma)
library(expsmooth)
library(nortest)
library(faraway)

### Ánálisis  ##################################################################################
## Identificación del modelo ####################################################################

# Primero, importamos la serie de tiempo
datos <- read.csv("Henry_Hub_Natural_Gas_Spot_Price.csv", header = TRUE)
datos <- datos[rev(rownames(datos)),]
NGSP <- ts(datos[,2], start=1997, freq=12)

# Usar estos comandos para verificar el directorio actual y cambiarlo, de ser necesario
# getwd()
# setwd("C:/Users/orteg/OneDrive/Documents/GitHub/AnalisisSeriesTiempo")

# Ahora graficaremos nuestra serie:
autoplot(NGSP, xlab = "Año", ylab = "Precio") 

# A simple vista, parece que tenemos un proceso no estacionario y sin ninguna tendencia estacional.
# Para veirificarlo, usamos los siguientes comandos, respectivamente:
adf.test(NGSP)
pp.test(NGSP)
?ndiffs
ggseasonplot(NGSP, xlab = "Mes",
             main = "",
             title = "asdfasdf", 
             season.labels = NULL,
             year.labels = FALSE,
             year.labels.left = FALSE,
             col = NULL,
             continuous = FALSE,
             polar = FALSE,
)

# Ambas pruebas, la de Dickey-Fuller aumenatada y la prueba de Phillips-Perron nos arrojan 
# valores p de 0.1607 y 0.1212, respectivamente, lo que nos indica que lidiamos con una serie
# no estacional.

# Para verificar si es necesario aplicar diferencias para convertir el proceso en uno estacionario
# y uno no estacional, revisaremos los resultados de los comandos siguientes, respectivamente:
ndiffs(NGSP)
nsdiffs(NGSP)

# Despues, observaremos cómo se comportan su función de autocorrelación (FAC) y función
# de autocorrelación parcial (FACP).
par(mfrow=c(2,1))
acf(NGSP, lag.max = 1000, main="FAC", ylab="FAC")
pacf(NGSP, lag.max = 1000, main="FACP", ylab="FACP")
par(mfrow=c(1,1))
# Al observar la FAC y la FACP podemos notar que el compartamiento se asemeja a un AR(1).
# Sin embargo, el proceso debe ser diferenciado (y posiblemente transformado) para ser estacionario,
# por lo que un posible modelo sea un ARIMA(1,1,0)

# Ahora, revisaremos lo que la función auto.arima nos devuelve
auto.arima(NGSP)
# La propuesta de modelo con la función auto.arima nos arroja un ARIMA(0,1,0), o un (I(1)). Es 
# decir, parece que tenemos un proceso de caminata aleatoria, la cual sabemos que no es estacionaria
# sino hasta después de aplicarle una primera diferencia.

# Después de este pequeño y breve análisis, tenemos dos posibles modelos
# ARIMA(1,1,0)
# vs
# ARIMA(0,1,0)

# Tener estos dos modelos similares nos lleva a pensar que el parámetro phi del ARIMA(1,1,0)
# será muy pequeño, por lo que tendremos especial cuidado en el supuesto de parsimonía para
# este parámetro.

# Estabilización de varianza ###################################################################
# Con el propósito de volver estacionaria la serie, buscaremos una transformación al proceso
# original para ver si ésta puede volverla una serie estacionaria.






# Probaremos con la metodología descrita en el capítulo 4 del libro "Análisis estadístico de 
# Series de Tiempo Económicas", de Victor Guerrero:
H <- 5
N <- length(NGSP)
n <- 9
R <- (N-n)/H

lambdas <- c(-9,-8.5,-8,-7.5,-7,-6.5,-6,5.5)
lambdas <- seq(from=-8.5,to=8.5,by =0.1)

breaks <- c()
breaks[1] <- 0
for(i in 1:H){
  breaks[i+1] <- R*i
}

temp <- as.vector(t(NGSP))
temp <- head(temp,length(temp)-n)

matrix <- matrix(rep(0), nrow = H, ncol = length(lambdas))
S <- c()
Z <- c()
for(i in 1:H){
  a <- breaks[i]+1
  b <- breaks[i+1]
  temp_2 <- temp[a:b]
  Z[i] <- sum(temp_2)/R
  S[i] <- sqrt((sum(temp_2-Z[i])^(2))/(R-1))
}


for(i in 1:H){
  for(j in 1:length(lambdas)){
    matrix[i,j] <- (S[i])/(Z[i]^(1-lambdas[j]))
  }
}

CC <- c()
M_lambda <- c()
de_lambda <- c()
for(i in 1:length(lambdas)){
  M_lambda[i] <- sum(matrix[,i])/H
  de_lambda[i] <- sqrt(sum(((matrix[,i])-(M_lambda[i]))^2)/(H-1))
  CC[i] <- de_lambda[i]/M_lambda[i]
}

CC <- rbind(lambdas,CC)
lambda <- CC[1,which.min(CC[2,])]

lambda

# De acuerdo a este análisis, parece que la transformación es de la forma f(x) = x^2.
# Graficaremos esta transformación para obtener una vista previa del proceso:
NGSPT <- NGSP^lambda
autoplot(NGSPT, xlab = "Tiempo", ylab = "Precio")
acf(NGSPT,lag.max = 100, xlab = "Lag", ylab = "FAC")
# La gráfica parece indicarnos que el proceso sigue siendo no estacionario. 

# Verificaremos con algunas prueba para confirmar estacionariedad.
adf.test(NGSPT)
pp.test(NGSPT)
autoplot(NGSPT)
auto.arima(NGSP)

# Revisaremos su FAC y FACP:
acf(NGSPT)
pacf(NGSPT)

# A pesar que la prueba de Dickey-Fuller aumentada y la prueba de Phillips-Perron nos arroja 
# un valor-p mucho más pequeño (con el que podríamos rechazar la no-estacionariedad), esta 
# función parece incrementar los picos presentados en 2002, 2005 y 2008.
# De igual manera, su FAC y FACP nos muestran que una primera diferencia sigue siendo
# necesaria, por lo que tendremos en cuenta esta transformación mas no la usaremos.

# La anterior transformación parece no ser la más adecuada para estabilizar la varianza, pues
# agranda los picos presentes en los años 2002, 2006 y 2008.

# Probaremos ahora con la transformación de Box-Cox, una de las transformaciones más conocidas:
lambda <- BoxCox.lambda(NGSP)
lambda_bc <- lambda
NGSP_BC <- BoxCox(NGSP, lambda)
autoplot(NGSP_BC, xlab="Tiempo", ylab="Precio")
# Parece que Box-Cox es una transformación util, pues se puede notar que la varianza es similar
# en diferentes instancias del tiempo. 
adf.test(NGSP_BC)
pp.test(NGSP_BC)
# Revisaremos su FAC y FACP:
acf(NGSP_BC, lag.max = 1000, main="FAC", ylab="FAC")
pacf(NGSP_BC, lag.max = 1000, main="FACP", ylab="FACP")
# Como habíamos comentado anteriormente, se asemejan al de un AR(1,1,0) o un I(1)
# Sin embargo, el proceso sigue siendo no estacionario.

# Al principio, utilizamos la serie original sin ninguna transformación estabilziadora de varianza.
# Sin embargo, el supuesto de independencia de residuos no se cumple, por lo que usaremos Box-Cox.
backup <- NGSP
NGSP <- NGSP_BC 

# Ahora utilizaremos el siguiente comando para ver cuántas
# diferencias son necesarias para volver el proceso a uno estacionario.
ndiffs(NGSP)

# Aplicamos entonces una primera diferencia a la serie original
DNGSP <- diff(NGSP)
autoplot(DNGSP)
acf(diff(NGSP), lag.max = 1000, main="FAC", ylab="FAC")
pacf(diff(NGSP), lag.max = 1000, main="FACP", ylab="FACP")

# Parece que la serie es ahora estacionaria. Para confirmar...
adf.test(DNGSP)
pp.test(DNGSP)
# ¿Punto de equilibrio?
mean(diff(NGSP_BC))
mean(diff(NGSP_BC)) + sqrt(var(diff(NGSP_BC)))
mean(diff(NGSP_BC)) - sqrt(var(diff(NGSP_BC)))
# En efecto, ambas pruebas arrojan valores p menores a 0.01
# Finalmente, revisamos su FAC y FACP
acf(diff(NGSP_BC), lag.max = 1000, main="FAC", ylab="FAC")
pacf(diff(NGSP_BC), lag.max = 1000, main="FACP", ylab="FACP")

# Ahora, la serie ya es estacionaria, pero ¿qué pasa si aplicamos una segunda diferencia?
acf(diff(diff(NGSP)), lag.max = 1000, main="FAC", ylab="FAC")
pacf(diff(diff(NGSP)), lag.max = 1000, main="FACP", ylab="FACP")
acf(diff(diff(NGSP_BC)), lag.max = 1000, main="FAC", ylab="FAC")
pacf(diff(diff(NGSP_BC)), lag.max = 1000, main="FACP", ylab="FACP")
# Por las graficas de arriba, percibimos un MA(0,2,1)
# Hay que tener cuidado con la sobrediferenciación, pues nuestro modelo
# podría perder poder predictivo


################################################################################################
# Los modelos que proponemos al observar su FAC y FACP son
### ARIMA(0,1,0) ### Pues se asemeja a la FAC y FACP de una caminata aleatoria
### ARIMA(1,1,0) ### Pues es p
### ARIMA(5,1,9) ### Por el modelo propuesto por Hosseinipoor



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
## Estimación de parámetros  ###################################################################
# ARIMA(0,1,0) #################################################################################
# model <- Arima(NGSP_BC, order=c(0,1,0), include.drift = TRUE)
# Verificando el valor de la deriva, podemos ver que es bastante pequeña, por lo que decidimos
# excluirla del modelo
model <- Arima(NGSP_BC, order=c(0,1,0), include.drift = FALSE) 

# Residuales
residuals <- residuals(model)
checkresiduals(model)


# Supuesto 1: Media cero #######################################################################
# Debemos verificar que el valor absoluto del cociente sea menor que dos para decir que no hay 
# evidencia de quela media del proceso sea diferente de 0. 
N <- length(residuals)
p <- 0
d <- 1
q <- 0

mean <- sum(residuals[-1])/(N-d-p)
std <- sqrt(sum((residuals[-1] - mean)^2)/(N-d-p-q))

cociente <- (sqrt(N-d-p)) * (mean/std)
abs(cociente)
# Como el valor absoluto del conciente es menor que 2, entonces podemos decir que no hay evidencia
# suficiente para afirmar que la media del proceso es distinta de cero.


# Supuesto 2: Varianza constante ###############################################################
# Observamos de manera visual si la varianza parece ser constante o no
checkresiduals(model)
# De forma visual, parece ser que la varianza parece ser constante, pues no se presenta mucha
# variación en estos a través de los datos históricos.


# Supuesto 3: Residuos independientes ##########################################################
# Aplicamos la prueba de Ljung-Box
# Un artículo interesante sobre cómo elegir el número de lags
# https://robjhyndman.com/hyndsight/ljung-box-test/
checkresiduals(model)
acf(residuals)
Box.test(residuals, type = c("Ljung-Box"), lag = 10)
Box.test(residuals, type = c("Box-Pierce"), lag = 10)
# Como podemos ver, ambas pruebas no rechazan el supuesto de independencia con la serie,
# transformada, caso contrario a la serie normal.


# Supuesto 4: Normalidad #######################################################################
# Verificar que aprox. el 95% de las observaciones se encuentren dentro del intervalo que se extiende 2 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
std <- sqrt(var(residuals))
left <- -2*std
right <- 2*std

temp <- length(residuals[residuals>right]) + length(residuals[residuals<left])
temp/length(residuals)
# Alrededor del 6% de los residuales están fuera del intervalo de dos desviaciones estándar, pero
# parece que este supuesto lo pasa, pues esperábamos un 5%.
# De igual forma, revisaremos si los residuales se distribuyen de forma simétrica. 
breaks <- seq(from = floor(max(residuals))-1, to = ceiling(max(residuals)), by = 0.1)
hist(residuals, breaks = breaks)
# Es simétrico, sin duda.

# Ahora, aplicaremos las pruebas de Lilliefors ni Shapiro-Wilk
shapiro.test(residuals)
lillie.test(residuals)
# No pasa las pruebas, pero es suficiente con las dos anteriores

# Además, al análizar el gráfico qq, algunas observaciones no caen dentro de la línea recta.
qqnorm(residuals)
qqline(residuals)


# Supuesto 5: No observaciones aberrantes ######################################################
# Prácticamente todas las observaciones deberían estar dentro del intervalo que se extiende 3 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
std <- sqrt(var(residuals))
left <- -3*std
right <- 3*std

temp <- length(residuals[residuals>right]) + length(residuals[residuals<left])
temp/length(residuals)
# Prácticamente todos los residuos están dentro del rango de tres desviaciones estándar.


# Supuesto 6: Parsimonía #######################################################################
# Ver con un 95% de confianza que todos los parámetros sean diferentes de 0.
# En este caso no hay parámetros ya que el modelo es I(1).


# Supuesto 7: Modelo admisible #################################################################
# Verificar que los parámetros se encuentren dentro de las regiones admisibles correspondientes.
# En este caso no hay parámetros ya que el modelo es I(1).


# Supuesto 8: Modelo estable ###################################################################
# Calculamos las correlaciones entre pares para ver que sean bajas.
# En este caso no hay correlaciones ya que no hay parámetros, debido a que el modelo es I(1).


## Pronósticos #################################################################################
forecast <- rwf(backup, drift=FALSE, lambda=lambda_bc, h=12, level=80,
           biasadj=TRUE)
autoplot(backup, xlab="Tiempo", ylab="Precios") +
  autolayer(forecast)
# Estos pronósticos parecen tener sentido, puedes nuesto modelo es un I(1).
# La demostración se deja al lector.



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
# ARIMA(5,1,9) #################################################################################
model <- Arima(NGSP_BC, order=c(5,1,9))
#?Arima
# Residuales
residuals <- residuals(model)
checkresiduals(model)


# Supuesto 1: Media cero #######################################################################
# Debemos verificar que el valor absoluto del cociente sea menor que dos para decir que no hay 
# evidencia de quela media del proceso sea diferente de 0. 
N <- length(residuals)
p <- 5
d <- 1
q <- 9

mean <- sum(residuals[-1])/(N-d-p)
std <- sqrt(sum((residuals[-1] - mean)^2)/(N-d-p-q))

cociente <- (sqrt(N-d-p)) * (mean/std)
abs(cociente)
# Como el valor absoluto del conciente es menor que 2, entonces podemos decir que no hay evidencia
# suficiente para afirmar que la media del proceso es distinta de cero.


# Supuesto 2: Varianza constante ###############################################################
# Observamos de manera visual si la varianza parece ser constante o no
checkresiduals(model)
# De forma visual, parece ser que la varianza parece ser constante, pues no se presenta mucha
# variación en estos a través de los datos históricos.


# Supuesto 3: Residuos independientes ##########################################################
# Aplicamos la prueba de Ljung-Box
checkresiduals(model)
acf(residuals)
Box.test(residuals, type = c("Ljung-Box"), lag = 10)
Box.test(residuals, type = c("Box-Pierce"), lag = 10)
# Como podemos ver, ambas pruebas no rechazan el supuesto de independencia con la serie,
# transformada, caso contrario a la serie normal.


# Supuesto 4: Normalidad #######################################################################
# Verificar que aprox. el 95% de las observaciones se encuentren dentro del intervalo que se extiende 2 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
std <- sqrt(var(residuals))
left <- -2*std
right <- 2*std

temp <- length(residuals[residuals>right]) + length(residuals[residuals<left])
temp/length(residuals)
# Alrededor del 6% de los residuales están fuera del intervalo de dos desviaciones estándar, pero
# parece que este supuesto lo pasa, pues esperábamos un 5%.
# De igual forma, revisaremos si los residuales se distribuyen de forma simétrica. 
breaks <- seq(from = floor(max(residuals))-1, to = ceiling(max(residuals)), by = 0.1)
hist(residuals, breaks = breaks)
# Es simétrico, sin duda.

# Ahora, aplicaremos las pruebas de Lilliefors ni Shapiro-Wilk
shapiro.test(residuals)
lillie.test(residuals)
# Estas pruebas sí las pasa, lo cual es genial.

# Además, al análizar el gráfico qq, algunas observaciones no caen dentro de la línea recta.
qqnorm(residuals)
qqline(residuals)


# Supuesto 5: No observaciones aberrantes ######################################################
# Prácticamente todas las observaciones deberían estar dentro del intervalo que se extiende 3 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
std <- sqrt(var(residuals))
left <- -3*std
right <- 3*std

temp <- length(residuals[residuals>right]) + length(residuals[residuals<left])
temp/length(residuals)
# Prácticamente todos los residuos están dentro del rango de tres desviaciones estándar.


# Supuesto 6: Parsimonía #######################################################################
# Ver con un 95% de confianza que todos los parámetros sean diferentes de 0.
mean <- model$coef
std <- c()
for(i in 1:length(mean)){
  std[i] <- sqrt(model$var.coef[i,i])
}
left <- mean - (2)*std
right <- mean + (2)*std
left*right
# De acuerdo a este pequeño análisis, es posible que los parámetros 1, 6 y 7 sean 0, debido a 
# que sus intervalos de confianza contienen al 0. Sin embargo, si reducimos un poco el nivel
# de significancia, este supuesto pueda cumplirse.


# Supuesto 7: Modelo admisible #################################################################
# Verificar que los parámetros se encuentren dentro de las regiones admisibles correspondientes.
# Si algo falta, es esta parte. Necesitaríamos calcular las determinantes de las matrices y usar
# las fórmulas de Jules-Walker.
# Otra forma de verificar este supuesto es simplemente confirmar estacionariedad e invertibilidad
ar <- mean[1:5]
ma <- mean[6:14]

install.packages("portes")
library(portes)
InvertQ(ar)
InvertQ(ma)
# Parece que el proceso no es invertible

# Como podemos observar, el modelo es invertible y estacionario, como ya habíamos verificado.
# Podemos decir que este supuesto se cumple.


# Supuesto 8: Modelo estable ###################################################################
# Calculamos las correlaciones entre pares para ver que sean bajas.
# Matriz de covarianzas
model$var.coef

model$var.coef[1,1]

std <- c()
for(i in 1:14){
  std[i] <- sqrt(model$var.coef[i,i])
}

temp <- model$var.coef
for(i in 1:14){
  for(j in 1:14){
    temp[i,j] <- temp[i,j]/(std[i]*std[j]) 
  }
}

(temp > 0.3)*temp
write.table((temp > 0.3)*temp, "C:/Users/orteg/OneDrive/Documents/GitHub/AnalisisSeriesTiempo/correlaciones.csv", sep="\t")


## Pronósticos #################################################################################
pronostico <-forecast(model, lambda = lambda_bc, h = 12, biasadj = TRUE)
autoplot(pronostico, main="", xlab="Tiempo", ylab="Precios") 

autoplot(backup, xlab="Tiempo", ylab="Precios") +
  autolayer(pronostico)

# Estos pronósticos parecen tener sentido, puedes nuesto modelo es un I(1).
# La demostración se deja al lector.



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
# ARIMA(0,2,1) #################################################################################
model <- Arima(NGSP_BC, order=c(0,2,1)) 

# Residuales
residuals <- residuals(model)
checkresiduals(model)


# Supuesto 1: Media cero #######################################################################
# Debemos verificar que el valor absoluto del cociente sea menor que dos para decir que no hay 
# evidencia de quela media del proceso sea diferente de 0. 
N <- length(residuals)
p <- 5
d <- 1
q <- 9

mean <- sum(residuals[-1])/(N-d-p)
std <- sqrt(sum((residuals[-1] - mean)^2)/(N-d-p-q))

cociente <- (sqrt(N-d-p)) * (mean/std)
abs(cociente)
# Como el valor absoluto del conciente es menor que 2, entonces podemos decir que no hay evidencia
# suficiente para afirmar que la media del proceso es distinta de cero.


# Supuesto 2: Varianza constante ###############################################################
# Observamos de manera visual si la varianza parece ser constante o no
checkresiduals(model)
# De forma visual, parece ser que la varianza parece ser constante, pues no se presenta mucha
# variación en estos a través de los datos históricos.


# Supuesto 3: Residuos independientes ##########################################################
# Aplicamos la prueba de Ljung-Box
checkresiduals(model)
acf(residuals)
Box.test(residuals, type = c("Ljung-Box"), lag = 10)
Box.test(residuals, type = c("Box-Pierce"), lag = 10)
# Como podemos ver, ambas pruebas no rechazan el supuesto de independencia con la serie,
# transformada, caso contrario a la serie normal.


# Supuesto 4: Normalidad #######################################################################
# Verificar que aprox. el 95% de las observaciones se encuentren dentro del intervalo que se extiende 2 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
std <- sqrt(var(residuals))
left <- -2*std
right <- 2*std

temp <- length(residuals[residuals>right]) + length(residuals[residuals<left])
temp/length(residuals)
# Alrededor del 6% de los residuales están fuera del intervalo de dos desviaciones estándar, pero
# parece que este supuesto lo pasa, pues esperábamos un 5%.
# De igual forma, revisaremos si los residuales se distribuyen de forma simétrica. 
breaks <- seq(from = floor(max(residuals))-1, to = ceiling(max(residuals)), by = 0.1)
hist(residuals, breaks = breaks)
# Es simétrico, sin duda.

# Ahora, aplicaremos las pruebas de Lilliefors ni Shapiro-Wilk
shapiro.test(residuals)
lillie.test(residuals)
# Estas pruebas sí las pasa, lo cual es genial.

# Además, al análizar el gráfico qq, algunas observaciones no caen dentro de la línea recta.
qqnorm(residuals)
qqline(residuals)


# Supuesto 5: No observaciones aberrantes ######################################################
# Prácticamente todas las observaciones deberían estar dentro del intervalo que se extiende 3 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
std <- sqrt(var(residuals))
left <- -3*std
right <- 3*std

temp <- length(residuals[residuals>right]) + length(residuals[residuals<left])
temp/length(residuals)
# Prácticamente todos los residuos están dentro del rango de tres desviaciones estándar.


# Supuesto 6: Parsimonía #######################################################################
# Ver con un 95% de confianza que todos los parámetros sean diferentes de 0.
mean <- model$coef
std <- c()
for(i in 1:length(mean)){
  std[i] <- sqrt(model$var.coef[i,i])
}
left <- mean - (2)*std
right <- mean + (2)*std
left*right
# De acuerdo a este pequeño análisis, es posible que los parámetros 1, 6 y 7 sean 0, debido a 
# que sus intervalos de confianza contienen al 0. Sin embargo, si reducimos un poco el nivel
# de significancia, este supuesto pueda cumplirse.
autoplot(diff(NGSP_BC), xlab="Tiempo", ylab="Precio")
# Supuesto 7: Modelo admisible #################################################################
# Verificar que los parámetros se encuentren dentro de las regiones admisibles correspondientes.
# El parámetro se encuentra dentro de la región admisible, pues su valor absoluto es menor a 1


# Supuesto 8: Modelo estable ###################################################################
# En este caso, como sólo hay un parámetro


## Pronósticos #################################################################################
pronostico <-forecast(model, h = 12, biasadj = TRUE)
autoplot(pronostico, main="", xlab="Tiempo", ylab="Precios")
# Estos pronósticos parecen tener sentido, puedes nuesto modelo es un I(1).
# La demostración se deja al lector.



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
# ARIMA(1,1,0) #################################################################################
model <- Arima(NGSP, order=c(1,1,0)) 

# Residuales
residuals <- residuals(model)
checkresiduals(model)


# Supuesto 1: Media cero #######################################################################
# Debemos verificar que el valor absoluto del cociente sea menor que dos para decir que no hay 
# evidencia de quela media del proceso sea diferente de 0. 
N <- length(residuals)
p <- 1
d <- 1
q <- 0

mean <- sum(residuals[-1])/(N-d-p)
std <- sqrt(sum((residuals[-1] - mean)^2)/(N-d-p-q))

cociente <- (sqrt(N-d-p)) * (mean/std)
abs(cociente)
# Como el valor absoluto del conciente es menor que 2, entonces podemos decir que no hay evidencia
# suficiente para afirmar que la media del proceso es distinta de cero.


# Supuesto 2: Varianza constante ###############################################################
# Observamos de manera visual si la varianza parece ser constante o no
checkresiduals(model)
# De forma visual, parece ser que la varianza parece ser constante, pues no se presenta mucha
# variación en estos a través de los datos históricos.


# Supuesto 3: Residuos independientes ##########################################################
# Aplicamos la prueba de Ljung-Box
checkresiduals(model)
acf(residuals)
Box.test(residuals, type = c("Ljung-Box"), lag = 10)
Box.test(residuals, type = c("Box-Pierce"), lag = 10)
# Como podemos ver, ambas pruebas no rechazan el supuesto de independencia con la serie,
# transformada, caso contrario a la serie normal.


# Supuesto 4: Normalidad #######################################################################
# Verificar que aprox. el 95% de las observaciones se encuentren dentro del intervalo que se extiende 2 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
std <- sqrt(var(residuals))
left <- -2*std
right <- 2*std

temp <- length(residuals[residuals>right]) + length(residuals[residuals<left])
temp/length(residuals)
# Alrededor del 6% de los residuales están fuera del intervalo de dos desviaciones estándar, pero
# parece que este supuesto lo pasa, pues esperábamos un 5%.
# De igual forma, revisaremos si los residuales se distribuyen de forma simétrica. 
breaks <- seq(from = floor(max(residuals))-1, to = ceiling(max(residuals)), by = 0.1)
hist(residuals, breaks = breaks)
# Es simétrico, sin duda.

# Ahora, aplicaremos las pruebas de Lilliefors ni Shapiro-Wilk
shapiro.test(residuals)
lillie.test(residuals)
# No pasa las pruebas, pero es suficiente con las dos anteriores

# Además, al análizar el gráfico qq, algunas observaciones no caen dentro de la línea recta.
qqnorm(residuals)
qqline(residuals)


# Supuesto 5: No observaciones aberrantes ######################################################
# Prácticamente todas las observaciones deberían estar dentro del intervalo que se extiende 3 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
std <- sqrt(var(residuals))
left <- -3*std
right <- 3*std

temp <- length(residuals[residuals>right]) + length(residuals[residuals<left])
temp/length(residuals)
# Prácticamente todos los residuos están dentro del rango de tres desviaciones estándar.


# Supuesto 6: Parsimonía #######################################################################
# Ver con un 95% de confianza que todos los parámetros sean diferentes de 0.
mean <- model$coef
std <- sqrt(model$var.coef)
mean - (2)*std
mean + (2)*std
# Para bien o para mal, el supuesto de parsimonia no se cumple en este caso, lo que nos lleva
# a pensar que el modelo anterior (ARIMA(0,1,0)) era suficiente para describir el comportamiento
# de la serie.


# Supuesto 7: Modelo admisible #################################################################
# Verificar que los parámetros se encuentren dentro de las regiones admisibles correspondientes.
(abs(mean)<1)
# En este caso, el modelo es admisible sin duda.


# Supuesto 8: Modelo estable ###################################################################
# Calculamos las correlaciones entre pares para ver que sean bajas.
# En este caso no hay correlaciones ya que sólo hay un parámetro.


## Pronósticos #################################################################################
pronostico <-forecast(model, lambda = lambda_bc, h = 12, biasadj = TRUE)
autoplot(pronostico, main="", xlab="Tiempo", ylab="Precios") 

autoplot(backup, xlab="Tiempo", ylab="Precios") +
  autolayer(pronostico)

# Estos pronósticos parecen tener sentido, puedes nuesto modelo es un I(1).
# La demostración se deja al lector.