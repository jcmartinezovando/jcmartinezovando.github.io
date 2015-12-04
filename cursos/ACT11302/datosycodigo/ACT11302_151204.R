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
# install.packages("actuar")
# install.packages("fExtremes")
# install.packages("fitdistrplus")
# install.packages("ggplot")

##	Masked:
# install.packages("fGarch")
# install.packages("fTrading")
# install.packages("timeDate")
# install.packages("timeSeries")
# install.packages("fBasics")
# install.packages("survival")
# install.packages("splines")

library("actuar")
library("fExtremes")
library("fitdistrplus")

path.datos <- "C:/JCMO.Academia/@Cursos/2015-II_Calculo Actuarial III/_datos"

path.code <- "C:/JCMO.Academia/@Cursos/2015-II_Calculo Actuarial III/_codigo"

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

# 		D)	Gumbel
fit.gumbel <- gumbelFit(xdatos)
summary(fit.gumbel)

# 		E)	Generalized Extreme Value
fit.gev <- gevFit(xdatos)
summary(fit.gev)

# 		F)	Generalized Pareto 
fit.gpd <- gpdFit(xdatos)
summary(fit.gpd)

# 		Comparación (Gráfica)
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
qqcomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
cdfcomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)
ppcomp(list(fit.weibull, fit.lnorm, fit.gamma), legendtext = plot.legend)


#	Matriz
is.matrix(datos)
datos <- as.matrix(datos)
 
#	Data Frame
is.data.frame(datos)
datos <- as.data.frame(datos)


#	#       Graficacion
#	ggplot(datos1, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)

#	x <- datos$LossinDKM
#	hillPlot(x, start = 15, ci = 0.95, doplot = TRUE, plottype = c("alpha", "xi"), labels = TRUE)

#	--------------------------------------------------------------------------
#			Distribución Agregada de Severidades
#	--------------------------------------------------------------------------

#	p.ej. Supongamos que la mejor distribución para las severidades individuales es Pareto Generalizada
#	Así, el momento t = 0 necesitamos generar la distribución para S(1)
#			S(1)=sum_{i=1}^{N(1)} X_i
#

#	--------------------------------------------------------------------------
#	Caso 1.- N(1) fijo...

N_1 <- 500 		#	Número de siniestros para el tiempo t=1
				#	Pensemos que tan solo la cartera de seguros de auto en México es de 6.5 millones de pólizas
M <- 1000		#	Número de simulaciones para S(1)	
xdatos_1 <- matrix(NaN, M, N_1)
i <- 1
for(i in 1:N_1){
  xdatos_1[,i] <- as.matrix(gpdSim(model = list(xi = 0.4915575, mu = 0, beta = 7.0403588), n = M, seed = NULL))
  }
S_1 <- matrix(NaN, M, 1)
m <- 1
for(m in 1:M){
  S_1[m] <- sum(xdatos_1[m,])
  }

#	Más gráficas
par(mfrow = c(2, 2))
hist(S_1,round(M/10))
emdPlot(S_1)
qqparetoPlot(S_1)
msratioPlot(S_1)

#	Descripción
summary(S_1)  
S_1_VaR <- VaR(S_1, alpha = 0.05, type = "sample", tail = c("lower", "upper"))
S_1_CVaR <- CVaR(S_1, alpha = 0.05, type = "sample", tail = c("lower", "upper"))

#	--------------------------------------------------------------------------
#	Caso 2.- N(1) aleatorio...

#	Supongamos que N(1) ~ Po(\lambda = 500)

lambda <- 500
N_1_sim <- rpois(M, lambda) 
S_1_sim <- matrix(NaN, M, 1)
m <-1
for(m in 1:M){
  xdatos_aux <- as.matrix(gpdSim(model = list(xi = 0.4915575, mu = 0, beta = 7.0403588), n = N_1_sim[m], seed = NULL))
  if(m==1){
	xdatos_1_sim <- list(xdatos_aux)
	}else{
	xdatos_1_sim <- list(xdatos_1_sim, xdatos_aux)
	} 
  S_1_sim[m] <- sum(xdatos_aux)
  }

#	Más gráficas
par(mfrow = c(2, 2))
hist(S_1_sim,round(M/10))
emdPlot(S_1_sim)
qqparetoPlot(S_1_sim)
msratioPlot(S_1_sim)

#	Descripción
summary(S_1_sim)  
S_1_VaR <- VaR(S_1_sim, alpha = 0.05, type = "sample", tail = c("lower", "upper"))
S_1_CVaR <- CVaR(S_1_sim, alpha = 0.05, type = "sample", tail = c("lower", "upper"))

#	--------------------------------------------------------------------------
#	Caso 2.- N(1) aleatorio...

#	Supongamos que N(1) ~ Po(\lambda = 500)

lambda <- 500
N_1_sim <- rpois(M, lambda) 
S_1_sim <- matrix(NaN, M, 1)
m <-1
for(m in 1:M){
  xdatos_aux <- as.matrix(gpdSim(model = list(xi = 0.4915575, mu = 0, beta = 7.0403588), n = N_1_sim[m], seed = NULL))
  if(m==1){
	xdatos_1_sim <- list(xdatos_aux)
	}else{
	xdatos_1_sim <- list(xdatos_1_sim, xdatos_aux)
	} 
  S_1_sim[m] <- sum(xdatos_aux)
  }

#	Más gráficas
par(mfrow = c(2, 2))
hist(S_1_sim,round(M/10))
emdPlot(S_1_sim)
qqparetoPlot(S_1_sim)
msratioPlot(S_1_sim)

#	Descripción
summary(S_1_sim)  
S_1_VaR <- VaR(S_1_sim, alpha = 0.05, type = "sample", tail = c("lower", "upper"))
S_1_CVaR <- CVaR(S_1_sim, alpha = 0.05, type = "sample", tail = c("lower", "upper"))

#	--------------------------------------------------------------------------
#	Caso 3.- N(1) aleatorio, con coaseguro (la compañía paga la proporción \alpha del siniestro)

#	Supongamos que N(1) ~ Po(\lambda = 500)

lambda <- 500
N_1_sim_coa <- rpois(M, lambda) 
S_1_sim_coa <- matrix(NaN, M, 1)
m <-1
alpha_coa <- 0.9
for(m in 1:M){
  xdatos_aux_coa <- as.matrix(gpdSim(model = list(xi = 0.4915575, mu = 0, beta = 7.0403588), n = N_1_sim_coa[m], seed = NULL))
  if(m==1){
	xdatos_1_sim_coa <- list(xdatos_aux_coa)
	}else{
	xdatos_1_sim_coa <- list(xdatos_1_sim_coa, xdatos_aux_coa)
	} 
  S_1_sim_coa[m] <- sum(xdatos_aux_coa)
  }

#	Más gráficas
par(mfrow = c(2, 2))
hist(S_1_sim_coa,round(M/10))
emdPlot(S_1_sim_coa)
qqparetoPlot(S_1_sim_coa)
msratioPlot(S_1_sim_coa)

#	Descripción
summary(S_1_sim_coa)  
S_1_VaR <- VaR(S_1_sim_coa, alpha = 0.05, type = "sample", tail = c("lower", "upper"))
S_1_CVaR <- CVaR(S_1_sim_coa, alpha = 0.05, type = "sample", tail = c("lower", "upper"))

#
#	-- FIN: ACT11302_151204.R --