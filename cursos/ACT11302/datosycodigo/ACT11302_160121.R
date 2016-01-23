#
#       ACT-11302: Calculo Actuarial III
#       
#       Autor:	Juan Carlos Martinez Ovando
#       Email:  juan.martinez.ovando@itam.mx
#
#       Este codigo fue probado en R v.3.1.2
#

rm(list = ls())

##	Main:
# install.packages("ggplot")

##	Masked:
 install.packages("actuar")
 install.packages("fExtremes")
 install.packages("fitdistrplus")
 install.packages("fGarch")
 install.packages("fTrading")
 install.packages("timeDate")
 install.packages("timeSeries")
 install.packages("fBasics")
 install.packages("survival")
 install.packages("splines")

library("actuar")
library("fExtremes")
library("fitdistrplus")

path.datos <- "C:/JCMO.Academia/@Cursos/2016-I_Calculo Actuarial III/_data&code"

#       Leemos los datos de siniestros individuales 
datos <- read.csv(
			paste(path.datos,"/act11302_DanishInsuranceData.csv", sep = ""), 
			header = TRUE)
# head(datos)
# tail(datos)
# colnames(datos)
# rownames(datos)

data(danishClaims)
write.csv(danishClaims, file = paste(path.datos,"/act11302_danishClaims.csv", sep = ""), row.names = FALSE)
xdatos <- danishClaims[,2]

# 	Estadisticas descriptivas
summary(xdatos)

#	descdist(xdatos, boot = 3000)

#	--------------------------------------------------------------------------
#			Análisis Descriptivo para Severidades
#	--------------------------------------------------------------------------
# 		Comparación (Gráfica)
plotdist(xdatos, histo = TRUE, demp = TRUE)

#	Más gráficas
par(mfrow = c(2, 2))
emdPlot(xdatos)
qqparetoPlot(xdatos)
msratioPlot(xdatos)

#
#	-- FIN: ACT11302_160121.R --