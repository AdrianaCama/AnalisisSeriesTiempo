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
NGSP <- tsclean(NGSP)
#autoplot(NGSP) 

# Propuesta de modelo con auto.arima sale que es ARIMA(0,1,0) (I(1)), o sea caminata aleatoria.
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


# Con el fin de convertir el proceso en uno estacionario, aplicaremos una transformación 
# estabilizadora de varianza. Primero, veremos Box Cox:
expected <- mean(NGSP)
expected
residuals <- residuals(naive(NGSP))
autoplot(NGSP)
autoplot(residuals)
checkresiduals(NGSP)

## Box Cox
transformacion <- BoxCox.lambda(NGSP)
BC <- BoxCox(NGSP, transformacion)
residuals <- residuals(naive(BC))
#autoplot(BC)
autoplot(residuals)

## Logaritmo
transformacion_2 <- log(NGSP)
residuals <- residuals(naive(transformacion_2))
autoplot(residuals)



# De acuerdo con esta página https://online.stat.psu.edu/stat462/node/148/, existen
# varias pruebas para revisar si la varianza de los residuales es constante. Generalmente es 
# suficiente verificarlo de manera visual pero existen pruebas que confirman nuestras sospechas


# Graficamos por año con el objetivo de observar si hay estacionalidad
ggseasonplot(NGSP)

# Variación de la seasonal plot con coordenadas polares. Se puede observar que no es estacional.
ggseasonplot(NGSP, polar = TRUE)

# Es posible que el resultado que nos arroja nsdiffs señale que el proceso no es estacional
nsdiffs(NGSP)

# Esta función es para ver el número de veces que hay que diferenciar una serie para volverla estacionaria, como 
# Como es una caminata aleatoria, sale 1, pues hemos demostrado anteriormente que la primera de la caminata aleatoria
# es un proceso estacionario.
ndiffs(NGSP)

# Ahora, aplicamos la primera diferencia para convertirlo a un proceso estacionario
NGSP_dif <- diff(NGSP, differences=1)





#Estabilización de varianza con transformación de Box-Cox
lambda0<-BoxCox.lambda(NGSP)
BoxCoxNGSP<-BoxCox(NGSP,lambda0)
autoplot(BoxCoxNGSP)


#Gráfica de funciones de autocorrelaci?n y test para verificar estacionariedad
ggAcf(BoxCoxNGSP, lag.max=100)
ggtsdisplay(BoxCoxNGSP, lag.max=100)

adf.test(BoxCoxNGSP) #Prueba de Dickey Fuller que no rechaza estacionariedad por lo que no habría que diferenciar







################# Análisis de residuos ##################
#Residuales
res_autoarima <- residuals(autoarima)
checkresiduals(autoarima)

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
# De forma visual, parece ser que la varianza no es constante desde 2002 hasta 2009, pues en esos años
# se presentan algunos picos que pueden afectar la varianza de los residuales. Para confirmar, ap



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