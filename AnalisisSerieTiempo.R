datos <- read.csv("Henry_Hub_Natural_Gas_Spot_Price.csv", header = TRUE)
datos <- datos[rev(rownames(datos)),]
NGSP <- ts(datos[,2], start=1997, freq=12)


# Propuesta de modelo con auto.arima sale que es ARIMA(0,1,0) (I(1)), o sea caminata aleatoria.
auto.arima(NGSP)

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
