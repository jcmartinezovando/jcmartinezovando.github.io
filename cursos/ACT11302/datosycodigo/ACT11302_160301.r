#
#       ACT-11302: Calculo Actuarial III
#       
#       Autor:	Juan Carlos Martinez Ovando
#       Email:  juan.martinez.ovando@itam.mx
#
#       Este codigo fue probado en R v.3.1.2
#
#		01 de marzo de 2016
#

rm(list = ls())

##	Main:
# install.packages("ggplot")

library("fExtremes")
library("fitdistrplus")
library("actuar")
library("fitdistrplus")
library("ggplot")

#	Windows
path <- "C:/JCMO.Academia/@Cursos/2015-II_Calculo Actuarial III/_data&code/"
path.plot <- "C:/JCMO.Academia/@Cursos/2016-I_Calculo Actuarial III/_notas/figuras/"

#	Linux
#path <- '/media/jcmo/ROCADATA/JCMO.Academia/@Cursos/2016-I_Calculo\ Actuarial\ III/JCMO.Academia/@Cursos/2015-II_Calculo\ Actuarial\ III/_data&code'


#       Leemos los datos de siniestros individuales 
datos <- read.csv(
			paste(path,"/act11302_DanishInsuranceData.csv", sep = ""), 
			header = TRUE)
head(datos)
tail(datos)
colnames(datos)
rownames(datos)

data(danishClaims)
write.csv(danishClaims, file = paste(path,"act11302_danishClaims.csv", sep = ""), row.names = FALSE)
xdatos <- danishClaims[,2]

# 	Estadisticas descriptivas
summary(xdatos)

#	descdist(xdatos, boot = 3000)

#	--------------------------------------------------------------------------
#			Análisis Descriptivo para Severidades
#	--------------------------------------------------------------------------
# 		Comparación (Gráfica)
pdf(paste(path.plot,"act11302_160301_fig01.pdf", sep = ""), width=7, height=5)
plotdist(xdatos, "norm", para=list(mean=mean(xdatos), sd=sd(xdatos)), histo = TRUE, demp = TRUE)
dev.off()

jpeg(paste(path.plot,"act11302_160301_fig01.jpg", sep = ""), width=700, height=500)
plotdist(xdatos, "norm", para=list(mean=mean(xdatos), sd=sd(xdatos)), histo = TRUE, demp = TRUE)
dev.off()

#	Más gráficas
pdf(paste(path.plot,"act11302_160301_fig02.pdf", sep = ""), width=7, height=5)
par(mfrow = c(2, 2))
hist(xdatos,100)
emdPlot(xdatos)
qqparetoPlot(xdatos)
msratioPlot(xdatos)
dev.off()

jpeg(paste(path.plot,"act11302_160301_fig02.jpg", sep = ""), width=700, height=500)
par(mfrow = c(2, 2))
hist(xdatos,100)
emdPlot(xdatos)
qqparetoPlot(xdatos)
msratioPlot(xdatos)
dev.off()

#	--------------------------------------------------------------------------
#			Estimación de Distribuciones para Severidades
#	--------------------------------------------------------------------------

# 		A)	Weibull
fit.weibull <- fitdist(xdatos, "weibull")
summary(fit.weibull)

# 		B)	Gamma
fit.gamma <- fitdist(xdatos, "gamma")
summary(fit.gamma)

# 		C)	Lognormal
fit.lnorm <- fitdist(xdatos, "lnorm")
summary(fit.lnorm)

# 		Comparación (Gráfica)
pdf(paste(path.plot,"act11302_160301_fig04.pdf", sep = ""), width=7, height=5)
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
qqcomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
cdfcomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
ppcomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
dev.off()

jpeg(paste(path.plot,"act11302_160301_fig04.jpg", sep = ""), width=700, height=500)
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
qqcomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
cdfcomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
ppcomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
dev.off()

# 		D)	Generalized Extreme Value
fit.gev <- gevFit(xdatos)
summary(fit.gev)
pdf(paste(path.plot,"act11302_160301_fig05.pdf", sep = ""), width=7, height=5)
plot(fit.gev)
dev.off()

jpeg(paste(path.plot,"act11302_160301_fig05.jpg", sep = ""), width=700, height=500)
summary(fit.gev)
dev.off()

# 		E)	Generalized Pareto 
fit.gpd <- gpdFit(xdatos)
summary(fit.gpd)
pdf(paste(path.plot,"act11302_160301_fig06.pdf", sep = ""), width=7, height=5)
plot(fit.gpd)
dev.off()

jpeg(paste(path.plot,"act11302_160301_fig06.jpg", sep = ""), width=700, height=500)
summary(fit.gpd)
dev.off()

#
#	-- FIN: ACT11302_150301.r --
