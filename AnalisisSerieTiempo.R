#install.packages("forecast")
library(forecast)
library(tseries)
library(fpp2)
library(ggplot2)
library(fma)
library(expsmooth)
library("nortest")

#package faraway
library(faraway)

################################################################################################
################################################################################################
################################################################################################
### Identificación del modelo ##################################################################
################################################################################################
################################################################################################
################################################################################################

# Dataset and Setup
# Primero, importamos la serie de tiempo
datos <- read.csv("Henry_Hub_Natural_Gas_Spot_Price.csv", header = TRUE)
datos <- datos[rev(rownames(datos)),]
NGSP <- ts(datos[,2], start=1997, freq=12)
# Con esto podemos ver la serie de tiempo
#autoplot(NGSP) 

# Segundo, quitamos outliers o missing values. CHECAR SI ESTÁ BIEN USAR ESTA FUNCIÓN xd
# NGSP <- tsclean(NGSP)

# Graficamos nuestra serie para obtener una vista previa
autoplot(NGSP) 

# Propuesta de modelo con auto.arima sale que es ARIMA(0,1,0) (I(1)), o sea caminata aleatoria.
# Eso sale sin usar el comando tsclean(NGSP). Usándolo sale ARIMA(1,1,0)
autoarima <- auto.arima(NGSP)

################################################################################################
# FAC, FACP y varianza
FAC <- acf(NGSP)
FACP <- pacf(NGSP)
VarNGSP<-var(NGSP)
VarNGSP
# Al observar la FAC y la FACP podemos notar que el compartamiento se asemeja a un AR(1). Por otra
# parte, la función auto.arima nos propone un modelo ARIMA(0,1,0). Analizaremos ambos.


# Dickey-Fuller Aumentado para probar estacionariedad. El proceso es no estacionario
adf.test(NGSP)

# Graficamos por año con el objetivo de observar si hay estacionalidad
ggseasonplot(NGSP)

# Variación de la seasonal plot con coordenadas polares. Se puede observar que no es estacional.
ggseasonplot(NGSP, polar = TRUE)

# Es posible que el resultado que nos arroja nsdiffs señale que el proceso no es estacional
nsdiffs(NGSP)


# Con el fin de convertir el proceso en uno estacionario, aplicaremos una transformación
# estabilizadora de varianza. Primero, veremos Box Cox:
expected <- mean(NGSP)
expected
residuals <- residuals(naive(NGSP))
autoplot(NGSP)
autoplot(residuals)

## Box Cox
lambda <- BoxCox.lambda(NGSP)
NGSP_BC <- BoxCox(NGSP, lambda)
# Después de aplicar la transformación estabilizadora de varianza, checamos si ya es estacionario
adf.test(NGSP_BC)
# Creo que está peor jaja :(

# Ahora probaremos la transformación logarítmica
## Logaritmo
NGSP_LOG <- log(NGSP)
adf.test(NGSP_LOG)
autoplot(NGSP_LOG)

# # Sin aplicar tsclean para quitar outliers, los valores p quedan 0.1607, 0.324 y .276 para la
# # serie normal, aplicando BoxCox y logaritmo
# # Aplicando tsclean quedan 0.2778, 0.3471 y 0.3443
#
# Probaremos otras transformaciones
## Square Root
NGSP_SR <- sqrt(NGSP)
adf.test(NGSP_SR)
autoplot(NGSP_SR)
# Sin tsclean 0.2204, con tsclean 0,3155
#
## Cube Root
NGSP_CR <- (NGSP)^(1/3)
adf.test(NGSP_CR)
autoplot(NGSP_CR)
# Sin tsclean 0.2392, con tsclean 0.3261

## One over
NGSP_OO <- 1/NGSP
adf.test(NGSP_OO)
autoplot(NGSP_OO)

# Parece que ninguna logra hacerla estacionaria, por lo que probaremos con la
# metodología que está en el Capítulo 4:
H <- 5
N <- length(NGSP)
n <- 9
R <- (N-n)/H
#lambdas <- c(-20,-10,-7,-5,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,7,10,20)
lambdas <- c(-9,-8.5,-8,-7.5,-7,-6.5,-6,5.5)
lambdas <- seq(from=-8.5,to=-7.5,by =0.05)

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
CC

which.min(CC[2,])

# Parece que la transformación debe ser elevar al cuadrado sin eliminar outliers ni 
# valores extremos
NGSP_SQ <- NGSP^2
adf.test(NGSP_SQ)
autoplot(NGSP_SQ)

# Ahora, si quitamos los valores extremos y outliers de la serie original, la mejor
# transformación es elevar la serie ~ a la -7 (-8.15)
NGSP_SR <- NGSP^(-8.15)
adf.test(NGSP_SR)
autoplot(NGSP_SR)

# Como podemos ver, la prueba de Dickey Fuller nos arroja un valor p de 0.05508, con lo que
# podríamos rechazar la no estacionariedad de la serie con esta transformación.
# Hay que recordar que esto es sin usar el comando tsclean


# Ahora utilizaremos el siguiente comando para ver cuántas
# diferencias son necesarias para volver el proceso a uno estacionario.
ndiffs(NGSP_SR)

# Aplicamos entonces una primera diferencia a la serie original
DNGSP <- diff(NGSP_SR)
autoplot(DNGSP)

#### ESTO ES UNA PRUEBA
# lambda <- BoxCox.lambda(DNGSP)
# NGSP_BC <- BoxCox(DNGSP, lambda)
# adf.test(NGSP_BC)
# autoplot(NGSP_BC)
# mean(DNGSP)

# Ahora veremos cómo se comportan la FAC y la FACP
DFAC <- acf(DNGSP)
DFACP <- pacf(DNGSP)

ggtsdisplay(DNGSP, lag.max=200)
#Aqui nos quedamos

# Se asemejan a lo que habíamos visto antes. Sabemos que puede ser un ARIMA(0,1,0)
# De hecho, con tsclean sale un ARIMA(1,1,0)
# pero intentaremos con un AR(1)

model <- Arima(DNGSP, order=c(0,2,1))
model <- Arima(DNGSP, order=c(0,1,0))
model <- Arima(DNGSP, order=c(1,1,0))
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
res_autoarima <- residuals(autoarima)
checkresiduals(autoarima)

### TEMPORAL
res_autoarima <- residuals(model)


#########################################################################################################
# Supuesto 1 (media cero)
#########################################################################################################
# Debemos verificar que el valor absoluto del cociente sea menor que dos para decir que no hay evidencia de que 
# la media del proceso sea diferente de 0. 
media <- mean(res_autoarima)
media
desv <- sqrt(var(res_autoarima))
desv

N <- length(res_autoarima)
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
checkresiduals(autoarima)
# De forma visual, parece ser que la varianza no es constante pues
# se presentan algunos picos que pueden afectar la varianza de los residuales. 
# Para confirmar, aplicaremos una prueba:


# Supuesto 3 (residuos independientes)
# Prueba de Ljung-Box
#Box.Ljung.Test(res_autoarima)


# Supuesto 4 (normalidad)
# Verificar que aprox. el 95% de las observaciones se encuentren dentro del intervalo que se extiende 2 
# desviaciones estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
LI2 <- -2*desv
LS2 <- 2*desv

qqnorm(res_autoarima)
qqline(res_autoarima)
checkresiduals(autoarima)
shapiro.test(res_autoarima)
lillie.test(x = res_autoarima)

# Supuesto 5 (no observaciones aberrantes)
# Prácticamente todas las observaciones deberían estar dentro del intervalo que se extiende 3 desviaciones
# estándar por arriba y por debajo de la media, la cual esperamos que sea 0.
LI3 <- -3*desv
LS3 <- 3*desv


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
autoarima_pronostico <-forecast(autoarima, h = 3) # Duda con h = 3
autoplot(autoarima_pronostico)