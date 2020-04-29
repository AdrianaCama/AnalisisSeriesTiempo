#install.packages("forecast")
library(forecast)
library(tseries)
library(fpp2)
library(ggplot2)
library(fma)
library(expsmooth)
library("nortest")
library(faraway)
getwd()

#setwd("C:/Users/orteg/OneDrive/Documents/GitHub/AnalisisSeriesTiempo")
################################################################################################
################################################################################################
################################################################################################
### Identificación del modelo ##################################################################
################################################################################################
################################################################################################
################################################################################################

# Dataset and Setup
# Primero, importamos la serie de tiempo
setwd("C:/Users/orteg/OneDrive/Documents/GitHub/AnalisisSeriesTiempo")


datos <- read.csv("serie_final.csv", header = TRUE)
#datos <- datos[rownames(datos),]
NGSP <- ts(datos[,2])
backup <- ts(datos[,2])

# Graficamos nuestra serie para obtener una vista previa
autoplot(NGSP) 

# Propuesta de modelo con auto.arima sale que es ARIMA(0,2,0)
autoarima <- auto.arima(NGSP)
autoarima

# FAC, FACP y varianza
FAC <- acf(NGSP)
FACP <- pacf(NGSP)

# Al observar la FAC y la FACP podemos notar que el compartamiento se asemeja a un AR(1)

# Dickey-Fuller Aumentado para probar estacionariedad. El proceso es no estacionario
adf.test(NGSP)
pp.test(NGSP)

#### Definitivamente no es estacionaria #####

# En este tipo de gráficas es conveniente usar una transformación logarítmica, pues parece
# existir una especie de crecimiento exponencial. Lo haremos en un momento.

# Como es un fenómeno reciente, no existen suficientes datos para afirmar estacionalidad
# ni nada por el estilo.


# Con el fin de convertir el proceso en uno estacionario, aplicaremos una transformación
# estabilizadora de varianza. Primero, veremos la transformación logarítmica
LNGSP <- log(NGSP)
autoplot(LNGSP)
adf.test(LNGSP)
##############################
# La transformación no es suficiente

##############################

# La serie original contiene 0, lo cual causa algunos problemas
# lambda <- BoxCox.lambda(NGSP)
# NGSP_BC <- BoxCox(NGSP, lambda)
# adf.test(NGSP_BC)
# 
# 
# expected <- mean(NGSP)
# expected
# residuals <- residuals(naive(NGSP))
# autoplot(NGSP)
# autoplot(residuals)
# 
# ## Box Cox
# 
# # Después de aplicar la transformación estabilizadora de varianza, checamos si ya es estacionario
# adf.test(NGSP_BC)
# pp.test(NGSP_BC)
# # Creo que está peor jaja :(
# 
# # Ahora probaremos la transformación logarítmica
# ## Logaritmo
# NGSP_LOG <- log(NGSP)
# adf.test(NGSP_LOG)
# autoplot(NGSP_LOG)
# 
# # # Sin aplicar tsclean para quitar outliers, los valores p quedan 0.1607, 0.324 y .276 para la
# # # serie normal, aplicando BoxCox y logaritmo
# # # Aplicando tsclean quedan 0.2778, 0.3471 y 0.3443
# #
# # Probaremos otras transformaciones
# ## Square Root
# NGSP_SR <- sqrt(NGSP)
# adf.test(NGSP_SR)
# autoplot(NGSP_SR)
# # Sin tsclean 0.2204, con tsclean 0,3155
# #
# ## Cube Root
# NGSP_CR <- (NGSP)^(1/3)
# adf.test(NGSP_CR)
# autoplot(NGSP_CR)
# # Sin tsclean 0.2392, con tsclean 0.3261
# 
# ## One over
# NGSP_OO <- 1/NGSP
# adf.test(NGSP_OO)
# autoplot(NGSP_OO)
# 
# # Parece que ninguna logra hacerla estacionaria, por lo que probaremos con la
# # metodología que está en el Capítulo 4:
# H <- 5
# N <- length(NGSP)
# n <- 6
# R <- (N-n)/H
# #lambdas <- c(-20,-10,-7,-5,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,7,10,20)
# lambdas <- c(-9,-8.5,-8,-7.5,-7,-6.5,-6,5.5)
# lambdas <- seq(from=-8.5,to=8.5,by =0.1)
# 
# breaks <- c()
# breaks[1] <- 0
# for(i in 1:H){
#   breaks[i+1] <- R*i
# }
# 
# temp <- as.vector(t(NGSP))
# temp <- head(temp,length(temp)-n)
# 
# matrix <- matrix(rep(0), nrow = H, ncol = length(lambdas))
# S <- c()
# Z <- c()
# for(i in 1:H){
#   a <- breaks[i]+1
#   b <- breaks[i+1]
#   temp_2 <- temp[a:b]
#   Z[i] <- sum(temp_2)/R
#   S[i] <- sqrt((sum(temp_2-Z[i])^(2))/(R-1))
# }
# 
# 
# for(i in 1:H){
#   for(j in 1:length(lambdas)){
#     matrix[i,j] <- (S[i])/(Z[i]^(1-lambdas[j]))
#   }
# }
# 
# CC <- c()
# M_lambda <- c()
# de_lambda <- c()
# for(i in 1:length(lambdas)){
#   M_lambda[i] <- sum(matrix[,i])/H
#   de_lambda[i] <- sqrt(sum(((matrix[,i])-(M_lambda[i]))^2)/(H-1))
#   CC[i] <- de_lambda[i]/M_lambda[i]
# }
# 
# CC <- rbind(lambdas,CC)
# CC
# 
# lambda <- CC[1,which.min(CC[2,])]
# lambda
# 
# # Parece que la transformación debe ser elevar al cuadrado sin eliminar outliers ni 
# # valores extremos
# NGSP <- NGSP^lambda
# adf.test(NGSP)
# autoplot(NGSP)
#autoplot(NGSP)

# Ahora, si quitamos los valores extremos y outliers de la serie original, la mejor
# transformación es elevar la serie ~ a la -7 (-8.15)
# NGSP_SR <- NGSP^(-8.15)
# adf.test(NGSP_SR)
# autoplot(NGSP_SR)

# Como podemos ver, la prueba de Dickey Fuller nos arroja un valor p de 0.05508, con lo que
# podríamos rechazar la no estacionariedad de la serie con esta transformación.
# Hay que recordar que esto es sin usar el comando tsclean


# De


# Vamo a probar esta transformación
# TNGSP <- sign(NGSP)*log1p(abs(NGSP))
# adf.test(TNGSP)
# autoplot(TNGSP)
# 
# NGSP[1]
# autoplot(NGSPmayor)
# 
# backup <- NGSP
# NGSP <- NGSP_BC


# Ahora utilizaremos el siguiente comando para ver cuántas
# diferencias son necesarias para volver el proceso a uno estacionario.
ndiffs(NGSP)

# Aplicamos entonces una primera diferencia a la serie original

# Esto es aplicando sólo una diferencia y la transformación logarítmica
D1NGSP <- diff(LNGSP, differences = 1)
autoplot(NGSP)
autoplot(D1NGSP)
adf.test(D1NGSP)


# Esto es no aplicando ninguna transformación
DNGSP <- diff(NGSP, differences = 2)
autoplot(NGSP)
autoplot(DNGSP)

adf.test(DNGSP)

# Esto es aplicando la transformación logarítmica
DLNGSP <- diff(LNGSP, differences = 2)
autoplot(NGSP)
autoplot(DLNGSP)

adf.test(DLNGSP)

###############################
# La serie se vuelve estacionaria cuando aplicamos dos
# diferencias de nivel y potencialmente una transformación logarítmica
###############################

# Ahora veremos cómo se comportan la FAC y la FACP

#ggtsdisplay(DNGSP, lag.max=200)

ggtsdisplay(DNGSP, lag.max = 200)
ggtsdisplay(LNGSP, lag.max = 200) # AR(1)?
ggtsdisplay(D1NGSP, lag.max=200)
ggtsdisplay(DNGSP, lag.max=200)
ggtsdisplay(DLNGSP, lag.max=200) #MA(1)

auto.arima(DLNGSP)
# DFAC <- acf(DNGSP)
# DFACP <- pacf(DNGSP)

# ggAcf(DNGSP, lag.max=200)
# ggtsdisplay(DNGSP, lag.max=200)
#Aqui nos quedamos

# Se asemejan a lo que habíamos visto antes. Sabemos que puede ser un ARIMA(0,1,0)
# De hecho, con tsclean sale un ARIMA(1,1,0)
# pero intentaremos con un AR(1)
# autoplot(NGSP)
# model <- Arima(NGSP, order=c(0,1,0)) # Prometedor
# 
# model <- Arima(NGSP, order=c(0,1,1)) #Nope
# model <- Arima(NGSP, order=c(0,2,1)) # Prometedor
# model <- Arima(NGSP, order=c(1,1,0)) # Prometedor
#


#####################################################################################
#####################################################################################
# Verificación de supuestos para ARIMA(0,2,1), obtenido con auto.arima
model <- Arima(LNGSP, order=c(0,2,1))
residuals <- residuals(model)
# checkresiduals(model)
# lillie.test(x = residuals)
# mean(residuals)
# sqrt(var())

# autoplot(forecast(model,h=10)
# )


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

# Podemos construir un intervalo de confianza para la media usando Bootstrap

# Como el valor absoluto del conciente es menor que 2, entonces podemos decir que no hay evidencia
# suficiente para afirmar que la media del proceso es distinta de cero.
#########################################################################################################
#PASS



#########################################################################################################
# Supuesto 2 (varianza constante)
#########################################################################################################
# Observamos de manera visual si la varianza parece ser constante o no
checkresiduals(model)
# De forma visual, parece ser que la varianza no es constante pues
# se presentan algunos picos que pueden afectar la varianza de los residuales. 
# Para confirmar, aplicaremos una prueba:
#########################################################################################################
# PASS (De momento)



#########################################################################################################
# Supuesto 3 (residuos independientes)
#########################################################################################################
# Prueba de Ljung-Box
checkresiduals(model)
# Los residuos son independientes
# Con un valor de signficancia del 95%
Box.test(residuals,type="Ljung")
# Parece que los residuos son independientes
#########################################################################################################
# PASS

#########################################################################################################
# Supuesto 4 (normalidad)
#########################################################################################################
# Verificar que aprox. el 95% de las observaciones se encuentren dentro del intervalo que se extiende 2 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
desv <- sqrt(var(residuals))
desv

LEFT <- -2*desv
RIGHT <- 2*desv

JIJI <- length(residuals[residuals>RIGHT]) + length(residuals[residuals<LEFT])
JIJI
4/length(residuals)

# La primera prueba la pasa, sólo el 6.66% de las obsevaciones sobrepasa
# las dos desviaciones estándar
breaks <- seq(from = -1, to = 1, by = 0.05)
hist(residuals, breaks = breaks)

# La distribución de los residuales parece ser simétrica. Hay algunos valores atípicos
# pero parece normal
#########################################################################################################
# PASS

# Estas pruebas no las pasa :(
qqnorm(residuals)
qqline(residuals)
checkresiduals(model)
shapiro.test(residuals)
lillie.test(x = residuals)




#########################################################################################################
# Supuesto 5 (no observaciones aberrantes)
#########################################################################################################
# Prácticamente todas las observaciones deberían estar dentro del intervalo que se extiende 3 desviaciones
# estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
desv <- sqrt(var(residuals))
desv

LEFT <- -3*desv
RIGHT <- 3*desv

JIJI <- length(residuals[residuals>RIGHT]) + length(residuals[residuals<LEFT])
JIJI
JIJI/length(residuals)

# El 3.33% de las observaciones (2) sobrepasan el intervalo de más menos tres std

checkresiduals(model)
autoplot(LNGSP)
autoplot(NGSP)

# Hay que checar esto, pero implícitamnete se asumió que no hay observaciones aberrantes
# Después de todo, se trata del número acumulado de casos confirmados de COVID en México


# plot(res_autoarima, ylim=c(-2*10^(-3),6*10^(-4)))
# plot(res_autoarima)
# abline(h=LS2, col="red")
# abline(h=LI2, col="red")
# abline(h=LI3, col="blue")
# abline(h=LS3, col="blue")

#########################################################################################################
# Supuesto 6 (parsimonía)
#########################################################################################################
# Ver con un 95% de confianza que todos los parámetros sean diferentes de 0.
meanie <- model$coef
se <- sqrt(model$var.coef)

meanie - (2)*se
meanie + (2)*se

# El cero no se encuentra dentro del intervalo
#########################################################################################################
# PASS


#########################################################################################################
# Supuesto 7 (modelo admisible)
#########################################################################################################
# Verificar que los parámetros se encuentren dentro de las regiones admisibles correspondientes.

# En este caso, sólo tenemos un parámetro que sería
meanie
# El cual definitivamente se encuentra dentro de las regiones admisibles
# De hecho, el intervalo de confianza nos asegura que está dentro de -.82, -.23 con n 95% de confianza

#########################################################################################################
# PASS


#########################################################################################################
# Supuesto 8 (modelo estable)
#########################################################################################################
# Calculamos las correlaciones entre pares para ver que sean bajas.
# En este caso no hay correlaciones ya que sólo hay un parámetro
#########################################################################################################
# PASS



####################### Pronósticos #######################
# auto.arima
fc <-forecast(model, h = 5) # Duda con h = 3
fc$mean<-exp(fc$mean)
fc$upper<-exp(fc$upper)
fc$lower<-exp(fc$lower)
fc$x<-exp(fc$x)
autoplot(fc)
