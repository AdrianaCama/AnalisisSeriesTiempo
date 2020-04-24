datos <- read.csv("Henry_Hub_Natural_Gas_Spot_Price.csv", header = TRUE)
datos <- datos[rev(rownames(datos)),]
NGSP <- ts(datos[,2], start=1997, freq=12)


# Propuesta de modelo con auto.arima sale que es ARIMA(0,1,0) (I(1)), o sea caminata aleatoria.
autoarima <- auto.arima(NGSP)


# FAC, FACP y varianza
FAC <- acf(NGSP)
FACP <- pacf(NGSP)
VarNGSP<-var(NGSP)
VarNGSP

# Dickey-Fuller Aumentado para probar estacionariedad. Sí lo es
adf.test(NGSP)

# Graficamos por año con el objetivo de observar si hay estacionalidad
ggseasonplot(NGSP)

# Esta función es para ver el número de veces que hay que diferenciar una serie para volverla estacionaria, como 
# ya es estacionaria, sale 0
nsdiffs(NGSP)


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


# Supuesto 1 (media cero)
# Debemos verificar que el valor absoluto del cociente sea menor que dos para decir que no hay evidencia de que 
# la media del proceso sea diferente de 0. 
media <- mean(autoarima)
media
desv <- sqrt(var(res_autoarima))
desv

N <- length(res_autoarima)
p <- 0
d <- 1
q <- 0
cociente <- (sqrt(N-d-p)) * (media/desv)
abs(cociente)


# Supuesto 2 (varianza constante)
# Observamos de manera visual si la varianza parece ser constante o no
checkresiduals(autoarima)


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