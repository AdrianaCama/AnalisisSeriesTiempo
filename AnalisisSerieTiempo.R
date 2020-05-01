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
library("nortest")
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
acf(NGSP)
pacf(NGSP)
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
pp.test(NGSP)

# A pesar que la prueba de Dickey-Fuller aumentada nos arroja un valor-p mucho más pequeño 
# (con el que podríamos rechazar la no-estacionariedad), la prueba de Phillips-Perron, una versión
# modificada de la prueba anterior y que proporcionada una conclusión más robusta acerca de la 
# posible estacionariedad, no cambia su valor-p.

# La anterior transformación parece no ser la más adecuada para estabilizar la varianza, pues
# agranda los picos presentes en los años 2002, 2006 y 2008.

# Probaremos ahora con la transformación de Box-Cox, una de las transformaciones más conocidas:
lambda <- BoxCox.lambda(NGSP)
NGSP_BC <- BoxCox(NGSP, lambda)
autoplot(NGSP_BC, xlab="Tiempo", ylab="Precio")
# Parece que Box-Cox es una transformación util, pues se puede notar que la varianza es similar
# en diferentes instancias del tiempo. 
adf.test(NGSP_BC)
pp.test(NGSP_BC)

# Revisaremos su FAC y FACP:
acf(NGSP)
pacf(NGSP_BC)
# Como habíamos comentado anteriormente, se asemejan al de un AR(1,1,0) o un I(1)

# Ahora compararemos la FAC y FACP de la serie transformada con Box Cox y la de un proceso
# de caminata aleatoria:
par(mfrow=c(2,2))
acf(NGSP_BC, main = "Precios de gas natural (Box-Cox)", ylab="FAC")
pacf(NGSP_BC, main = "Precios de gas natural (Box-Cox)", ylab="FACP")
random_walk <- arima.sim(model = list(order = c(0, 1, 0)), n = 23*12)
acf(random_walk, main = "Caminata aleatoria", ylab="FAC")
pacf(random_walk, main = "Caminata aleatoria", ylab="FACP")
par(mfrow=c(1,1))

# Realizamos un proceso similar con la FAC y FACP de la serie transformada con Box Cox 
# y la de un proceso ARIMA(1,1,0) con su parámetro menor a 0.1:
par(mfrow=c(2,2))
acf(NGSP_BC, main = "Precios de gas natural (Box-Cox)", ylab="FAC")
pacf(NGSP_BC, main = "Precios de gas natural (Box-Cox)", ylab="FACP")
ar <- arima.sim(list(order=c(1,1,0), ar=.1), n=500)
acf(ar, main = expression(paste("Proceso AR(1) con ",phi,"< 0.1")), ylab="FAC")
pacf(ar, main = expression(paste("Proceso AR(1) con ",phi,"< 0.1")), ylab="FACP")
par(mfrow=c(1,1))

### PENDIENTE
# Potencialmente podríamos incluir esto
par(mfrow=c(2,2))
acf(NGSP_BC, main = "Precios de gas natural (Box-Cox)", ylab="FAC")
pacf(NGSP_BC, main = "Precios de gas natural (Box-Cox)", ylab="FACP")
ar1<- arima.sim(list(order=c(1,1,2), ar=0.63371, ma=c(0.40338,0.08182)), n=500)
acf(ar1, main = expression(paste("Proceso AR(1) con ",phi,"< 0.1")), ylab="FAC")
pacf(ar1, main = expression(paste("Proceso AR(1) con ",phi,"< 0.1")), ylab="FACP")
par(mfrow=c(1,1))

# Después de aplicar la transformación estabilizadora de varianza, checamos si ya es estacionario
adf.test(NGSP_BC)
pp.test(NGSP_BC)
# Creo que está peor jaja :(

backup <- NGSP
NGSP <- NGSP_BC

# Ahora utilizaremos el siguiente comando para ver cuántas
# diferencias son necesarias para volver el proceso a uno estacionario.
ndiffs(NGSP)

# Aplicamos entonces una primera diferencia a la serie original
DNGSP <- diff(NGSP)
autoplot(DNGSP)

## Estimación de parámetros  ###################################################################
# Para un I(1)
model <- Arima(NGSP, order=c(0,1,0)) # Este pasó :D
# model <- Arima(NGSP, order=c(1,1,0)) # Este sólo no pasa parsimonía
# model$coef
# 
# 
# model <- Arima(NGSP, order=c(0,1,0), include.drift = TRUE)
# model <- Arima(NGSP, order=c(0,1,1)) #Nope
# model <- Arima(NGSP, order=c(0,2,1)) # Prometedor
# model <- Arima(NGSP, order=c(1,1,0)) # Prometedor

residuals <- residuals(model)
checkresiduals(model)

# # Esta función es para ver el número de veces que hay que diferenciar una serie para volverla estacionaria, como 
# # Como es una caminata aleatoria, sale 1, pues hemos demostrado anteriormente que la primera de la caminata aleatoria
# # es un proceso estacionario.
# ndiffs(NGSP)

# # Ahora, aplicamos la primera diferencia para convertirlo a un proceso estacionario
# NGSP_dif <- diff(NGSP, differences=1)


# #Estabilización de varianza con transformación de Box-Cox
# lambda0<-BoxCox.lambda(NGSP)
# BoxCoxNGSP<-BoxCox(NGSP,lambda0)
# autoplot(BoxCoxNGSP)


################# Análisis de residuos ##################
#Residuales
res_autoarima <- residuals(model)
checkresiduals(model)


#########################################################################################################
# Supuesto 1 (media cero)
#########################################################################################################
# Debemos verificar que el valor absoluto del cociente sea menor que dos para decir que no hay evidencia de que 
# la media del proceso sea diferente de 0. 
media <- mean(residuals)
media
desv <- sqrt(var(residuals))
desv

N <- length(residuals)
p <- 0
d <- 1
q <- 0
cociente <- (sqrt(N-d-p)) * (media/desv)
abs(cociente)
# Como el valor absoluto del conciente es menor que 2, entonces podemos decir que no hay evidencia
# suficiente para afirmar que la media del proceso es distinta de cero.


#########################################################################################################
# Supuesto 2 (varianza constante)
#########################################################################################################
# Observamos de manera visual si la varianza parece ser constante o no
checkresiduals(model)
# De forma visual, parece ser que la varianza no es constante pues
# se presentan algunos picos que pueden afectar la varianza de los residuales. 
# Para confirmar, aplicaremos una prueba:
# En el 2006 parece haber un problema (buscar que pasó) y ponerlo en el escrito
# Idea de mi amorcito preciosa <3



# Supuesto 3 (residuos independientes)
# Prueba de Ljung-Box
checkresiduals(model)
Box.test(residuals, type = c("Ljung-Box"))

#NULL HYPOTHESIS: INDEPENDENCE
# Parecer que a partir de un lag  de 9 se empieza a rechazar independencia
## PREGUNTA
# https://stats.stackexchange.com/questions/6455/how-many-lags-to-use-in-the-ljung-box-test-of-a-time-series
# Potencialmente usar prueba de correlación de atrasos para checar independencia
# CHECAR 
# PREGUNTA


# Los residuos son independientes

# Supuesto 4 (normalidad)
# Verificar que aprox. el 95% de las observaciones se encuentren dentro del intervalo que se extiende 2 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
LI2 <- -2*desv
LS2 <- 2*desv

desv <- sqrt(var(residuals))
desv

LEFT <- -2*desv
RIGHT <- 2*desv

JIJI <- length(residuals[residuals>RIGHT]) + length(residuals[residuals<LEFT])
JIJI
JIJI/length(residuals)
# Se esperaba 5% y es lo que pasó lol

breaks <- seq(from = -10, to = 10, by = 0.5)
hist(residuals, breaks = breaks)
# Sí se ve simétrico lol

qqnorm(residuals)
qqline(residuals)
checkresiduals(model)
shapiro.test(residuals)
lillie.test(x = residuals)


# Supuesto 5 (no observaciones aberrantes)
# Prácticamente todas las observaciones deberían estar dentro del intervalo que se extiende 3 desviaciones
# estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
LI3 <- -3*desv
LS3 <- 3*desv

desv <- sqrt(var(residuals))
desv

LEFT <- -3*desv
RIGHT <- 3*desv

JIJI <- length(residuals[residuals>RIGHT]) + length(residuals[residuals<LEFT])
JIJI
JIJI/length(residuals)
# 3% que era más o menos lo que esperábamos



plot(res_autoarima, ylim=c(-2*10^(-3),6*10^(-4)))
plot(res_autoarima)
abline(h=LS2, col="red")
abline(h=LI2, col="red")
abline(h=LI3, col="blue")
abline(h=LS3, col="blue")




# Supuesto 6 (parsimonía)
# Ver con un 95% de confianza que todos los parámetros sean diferentes de 0.
# En este caso no hay parámetros ya que el modelo es I(1).

# Supuesto 7 (modelo admisible)
# Verificar que los parámetros se encuentren dentro de las regiones admisibles correspondientes.
# En este caso no hay parámetros ya que el modelo es I(1).

# Supuesto 8 (modelo estable)
# Calculamos las correlaciones entre pares para ver que sean bajas.
# En este caso no hay correlaciones ya que no hay parámetros, debido a que el modelo es I(1).




####################### Pronósticos #######################
# auto.arima
autoarima_pronostico <-forecast(autoarima, h = 5) # Duda con h = 3
autoplot(autoarima_pronostico)

#https://stats.stackexchange.com/questions/333092/why-i-get-the-same-predict-value-in-arima-model
autoarima_pronostico

# Los pronósticos parecen estar bien 
# Te amo <3



































model <- Arima(NGSP, order=c(1,1,0)) # Este sólo no pasa parsimonía


#########################################################################################################
# Supuesto 1 (media cero)
#########################################################################################################
# Debemos verificar que el valor absoluto del cociente sea menor que dos para decir que no hay evidencia de que 
# la media del proceso sea diferente de 0. 
media <- mean(residuals)
media
desv <- sqrt(var(residuals))
desv

N <- length(residuals)
p <- 0
d <- 1
q <- 0
cociente <- (sqrt(N-d-p)) * (media/desv)
abs(cociente)

## ARIMA(1,1,0) pasa
# Como el valor absoluto del conciente es menor que 2, entonces podemos decir que no hay evidencia
# suficiente para afirmar que la media del proceso es distinta de cero.


#########################################################################################################
# Supuesto 2 (varianza constante)
#########################################################################################################
# Observamos de manera visual si la varianza parece ser constante o no
checkresiduals(model)
# De forma visual, parece ser que la varianza no es constante pues
# se presentan algunos picos que pueden afectar la varianza de los residuales. 
# Para confirmar, aplicaremos una prueba:
# En el 2006 parece haber un problema (buscar que pasó) y ponerlo en el escrito
# Idea de mi amorcito preciosa <3
# Checar de nuevo xd


# Supuesto 3 (residuos independientes)
# Prueba de Ljung-Box
checkresiduals(model)
Box.test(residuals, lag = 9, type = c("Ljung-Box"))

# Test de significancia para ver si autocorrelaciones son 
# significativamente distintas de 0?

# Parecer que a partir de un lag  de 9 se empieza a rechazar independencia
## PREGUNTA
# https://stats.stackexchange.com/questions/6455/how-many-lags-to-use-in-the-ljung-box-test-of-a-time-series
# Potencialmente usar prueba de correlación de atrasos para checar independencia
# CHECAR 
# PREGUNTA


# Los residuos son independientes

# Supuesto 4 (normalidad)
# Verificar que aprox. el 95% de las observaciones se encuentren dentro del intervalo que se extiende 2 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
LI2 <- -2*desv
LS2 <- 2*desv

desv <- sqrt(var(residuals))
desv

LEFT <- -2*desv
RIGHT <- 2*desv

JIJI <- length(residuals[residuals>RIGHT]) + length(residuals[residuals<LEFT])
JIJI
JIJI/length(residuals)
# Se esperaba 5% y es lo que pasó lol

breaks <- seq(from = -10, to = 10, by = 0.5)
hist(residuals, breaks = breaks)
# Sí se ve simétrico lol

qqnorm(residuals)
qqline(residuals)
checkresiduals(model)
shapiro.test(residuals)
lillie.test(x = residuals)



# Supuesto 5 (no observaciones aberrantes)
# Prácticamente todas las observaciones deberían estar dentro del intervalo que se extiende 3 desviaciones
# estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
LI3 <- -3*desv
LS3 <- 3*desv

desv <- sqrt(var(residuals))
desv

LEFT <- -3*desv
RIGHT <- 3*desv

JIJI <- length(residuals[residuals>RIGHT]) + length(residuals[residuals<LEFT])
JIJI
JIJI/length(residuals)
# 3% que era más o menos lo que esperábamos



plot(res_autoarima, ylim=c(-2*10^(-3),6*10^(-4)))
plot(res_autoarima)
abline(h=LS2, col="red")
abline(h=LI2, col="red")
abline(h=LI3, col="blue")
abline(h=LS3, col="blue")




# Supuesto 6 (parsimonía)
# Ver con un 95% de confianza que todos los parámetros sean diferentes de 0.
# En este caso no hay parámetros ya que el modelo es I(1).
meanie <- model$coef
se <- sqrt(model$var.coef)

meanie - (2)*se
meanie + (2)*se

# Parece ser que el supuesto de parsimonía no se cumple
# Parece no ser tan grave, tiene sentido que sea modelado por un ARI(1,1)
# Justificar bien



# Supuesto 7 (modelo admisible)
# Verificar que los parámetros se encuentren dentro de las regiones admisibles correspondientes.
# En este caso no hay parámetros ya que el modelo es I(1).
(abs(meanie)<1)


# Supuesto 8 (modelo estable)
# Calculamos las correlaciones entre pares para ver que sean bajas.

# Como sólo hay un parámetro, no hay autocorrelaciones entre parámetros



####################### Pronósticos #######################
# auto.arima
autoarima_pronostico <-forecast(autoarima, h = 5) # Duda con h = 3
autoplot(autoarima_pronostico)

#https://stats.stackexchange.com/questions/333092/why-i-get-the-same-predict-value-in-arima-model
autoarima_pronostico


asdfasdf <- rwf(NGSP, drift=FALSE, h=50, level=80, biasadj=TRUE)
#INVESTIGAR BIEN ESTO DE BIAS ADJ
autoplot(NGSP) + autolayer(asdfasdf, series="xd") +  guides(colour=guide_legend(title="Forecast"))

# Los pronósticos parecen estar bien 
# Te amo <3