#
#       ACT-11302: Cálculo Actuarial III
#       
#       Autor:	Juan Carlos Martinez Ovando
#       Email:		juan.martinez.ovando@itam.mx
#
#		Sesión:	3 de Mayo de 2016
#
#       Este código fue probado en R v.3.2.4
#

rm(list = ls())

#			Packages
# install.packages("fExtremes")
# install.packages("fitdistrplus")
# install.packages("ggplot2")

library("fExtremes")
library("fitdistrplus")
library("ggplot2")

#	Windows
path <- "C:/JCMO.Academia/@Mis.Cursos/2016-I_Calculo Actuarial III/_data&code/"
path.plot <- "C:/JCMO.Academia/@Mis.Cursos/2016-I_Calculo Actuarial III/_notas/figuras/"

#	Linux
#path <- '/media/jcmo/ROCADATA/JCMO.Academia/@Cursos/2016-I_Calculo\ Actuarial\ III/JCMO.Academia/@Cursos/2015-II_Calculo\ Actuarial\ III/_data&code/'
#path <- '/media/jcmo/ROCADATA/JCMO.Academia/@Cursos/2016-I_Calculo\ Actuarial\ III/JCMO.Academia/@Cursos/2015-II_Calculo\ Actuarial\ III/_notas/figuras/'


#       Leemos los datos de siniestros individuales 
datos <- read.csv(paste(path,"ACT11302_ClaimPred.csv", sep = ""), header = TRUE)

head(datos)
colnames(datos)

#	--------------------------------------------------------------------------
#			Análisis Descriptivo de Severidades
#	--------------------------------------------------------------------------

# 		Comparación (Gráfica)
datos.sev.2005 <- as.matrix(datos[which((datos$Model_Year=="2005") & (datos$Claim_Amount>0 ) ),"Claim_Amount"] )
datos.sev.2006 <- as.matrix(datos[which((datos$Model_Year=="2006") & (datos$Claim_Amount>0 ) ),"Claim_Amount"] )
datos.sev.2007 <- as.matrix(datos[which((datos$Model_Year=="2007") & (datos$Claim_Amount>0 ) ),"Claim_Amount"] )

#		Histogramas
sev.2005 <- hist(datos.sev.2005, 50, main="Individual Severity", xlab="Data")
sev.2006 <- hist(datos.sev.2006, 50) 
sev.2007 <- hist(datos.sev.2007, 50)

plot( sev.2005, col=rgb(0,0,1,1/4), xlim=c(0,4000), main="Individual Severity", xlab="Data")
plot( sev.2006, col=rgb(1,0,0,1/4), xlim=c(0,4000), add=T) 
plot( sev.2007, col=rgb(0,1,0,1/4), xlim=c(0,4000), add=T)

pdf(paste(path.plot,"act11302_160503_fig01.pdf", sep = ""), width=7, height=5)
plot( sev.2005, col=rgb(0,0,1,1/4), xlim=c(0,4000), main="Individual Severity", xlab="Data")
plot( sev.2006, col=rgb(1,0,0,1/4), xlim=c(0,4000), add=T) 
plot( sev.2007, col=rgb(0,1,0,1/4), xlim=c(0,4000), add=T)
dev.off()

jpeg(paste(path.plot,"act11302_160503_fig01.jpg", sep = ""), width=700, height=500)
plot( sev.2005, col=rgb(0,0,1,1/4), xlim=c(0,4000), main="Individual Severity", xlab="Data")
plot( sev.2006, col=rgb(1,0,0,1/4), xlim=c(0,4000), add=T) 
plot( sev.2007, col=rgb(0,1,0,1/4), xlim=c(0,4000), add=T)
dev.off()

#	--------------------------------------------------------------------------
#		Distribución de Severidades
#	--------------------------------------------------------------------------

# 		2005 & 2006 & 2007
fit.sev.gev <- gevFit(c(datos.sev.2005,datos.sev.2006,datos.sev.2007))
fit.sev.gpd <- gpdFit(c(datos.sev.2005,datos.sev.2006,datos.sev.2007))

par(mfrow = c(2, 2))
summary(fit.sev.gev)

pdf(paste(path.plot,"act11302_160503_fig03_gev.pdf", sep = ""), width=7, height=5)
par(mfrow = c(2, 2))
summary(fit.sev.gev)
dev.off()

jpeg(paste(path.plot,"act11302_160503_fig03_gev.jpg", sep = ""), width=700, height=500)
par(mfrow = c(2, 2))
summary(fit.sev.gev)
dev.off()

par(mfrow = c(2, 2))
summary(fit.sev.gpd)

pdf(paste(path.plot,"act11302_160503_fig03_gpd.pdf", sep = ""), width=7, height=5)
par(mfrow = c(2, 2))
summary(fit.sev.gpd)
dev.off()

jpeg(paste(path.plot,"act11302_160503_fig03_gpd.jpg", sep = ""), width=700, height=500)
par(mfrow = c(2, 2))
summary(fit.sev.gpd)
dev.off()

sev.lik <- as.data.frame(matrix(NA,nrow=2,ncol=2))
colnames(sev.lik) <- c("Model","log-Likelihood")
sev.lik[1,"Model"] <- "GEV"
sev.lik[2,"Model"] <- "GP"
sev.lik[1,"log-Likelihood"] <- 5602.74
sev.lik[2,"log-Likelihood"] <- 369.6645 

write.csv(sev.lik, file = paste(path.plot,"act11302_160503_tab01.csv", sep = ""), 
			sep = ",", row.names = TRUE, col.names = TRUE)

summ.fit.sev.gpd <- summary(fit.sev.gpd)

sev.dist <- as.data.frame(matrix(NA,nrow=2,ncol=4))
colnames(sev.dist) <- c("Estimator","xi","mu","beta")
sev.dist[1,"Estimator"] <- "Mean"
sev.dist[2,"Estimator"] <- "StdDev"

sev.dist[1,"xi"] <- 1.370929 
sev.dist[1,"mu"] <- 20.052697
sev.dist[1,"beta"] <- 30.410379 

sev.dist[2,"xi"] <- 0.05458572 
sev.dist[2,"mu"] <- 1.13667417
sev.dist[2,"beta"] <- 1.83501744 

write.csv(sev.dist, file = paste(path.plot,"act11302_160503_tab02.csv", sep = ""),  sep = ",", row.names = TRUE, col.names = TRUE)

#	--------------------------------------------------------------------------
#			Análisis Descriptivo de Frecuencias
#	--------------------------------------------------------------------------
datos.frec <- as.data.frame(matrix(NA, ncol=4,nrow=4))
colnames(datos.frec) <- c("Anio","Siniestro", "NoSiniestro","Frecuencia")
datos.frec[1,"Anio"] <- 2005
datos.frec[2,"Anio"] <- 2006
datos.frec[3,"Anio"] <- 2007
datos.frec[4,"Anio"] <- "Total"

datos.frec[1,"Siniestro"] <- nrow(as.matrix(datos[which((datos$Model_Year=="2005") & (datos$Claim_Amount>0 ) ),"Claim_Amount"] ))
datos.frec[2,"Siniestro"] <- nrow(as.matrix(datos[which((datos$Model_Year=="2006") & (datos$Claim_Amount>0 ) ),"Claim_Amount"] ))
datos.frec[3,"Siniestro"] <- nrow(as.matrix(datos[which((datos$Model_Year=="2007") & (datos$Claim_Amount>0 ) ),"Claim_Amount"] ))
datos.frec[4,"Siniestro"] <- sum(datos.frec[1:3,"Siniestro"])

datos.frec[1,"NoSiniestro"] <- nrow(as.matrix(datos[which((datos$Model_Year=="2005") & (datos$Claim_Amount==0 ) ),"Claim_Amount"] ))
datos.frec[2,"NoSiniestro"] <- nrow(as.matrix(datos[which((datos$Model_Year=="2006") & (datos$Claim_Amount==0 ) ),"Claim_Amount"] ))
datos.frec[3,"NoSiniestro"] <- nrow(as.matrix(datos[which((datos$Model_Year=="2007") & (datos$Claim_Amount==0 ) ),"Claim_Amount"] ))
datos.frec[4,"NoSiniestro"] <- sum(datos.frec[1:3,"NoSiniestro"])

datos.frec[1,"Frecuencia"] <- datos.frec[1,"Siniestro"] / (datos.frec[1,"Siniestro"]  + datos.frec[1,"NoSiniestro"])
datos.frec[2,"Frecuencia"] <- datos.frec[2,"Siniestro"] / (datos.frec[2,"Siniestro"]  + datos.frec[2,"NoSiniestro"])
datos.frec[3,"Frecuencia"] <- datos.frec[3,"Siniestro"] / (datos.frec[3,"Siniestro"]  + datos.frec[3,"NoSiniestro"])
datos.frec[4,"Frecuencia"] <- datos.frec[4,"Siniestro"] / (datos.frec[4,"Siniestro"]  + datos.frec[4,"NoSiniestro"])

write.csv(datos.frec, file = paste(path.plot,"act11302_160503_tab03.csv", sep = ""), 
			sep = ",", row.names = TRUE, col.names = TRUE)

plot(ts(datos.frec[1:3,c("Siniestro","NoSiniestro","Frecuencia")]),main="Trends",xlab="Time")

pdf(paste(path.plot,"act11302_160503_fig02.pdf", sep = ""), width=7, height=5)
plot(ts(datos.frec[1:3,c("Siniestro","NoSiniestro","Frecuencia")]),main="Trends",xlab="Time")
dev.off()

jpeg(paste(path.plot,"act11302_160503_fig02.jpg", sep = ""), width=700, height=500)
plot(ts(datos.frec[1:3,c("Siniestro","NoSiniestro","Frecuencia")]),main="Trends",xlab="Time")
dev.off()

#	--------------------------------------------------------------------------
#			Distribución del Monto Agregado de Siniestros
#	--------------------------------------------------------------------------
Msim <- 100000
S_dist <- as.data.frame(matrix(NA, ncol = 3, nrow = Msim))
S_dist <- as.matrix(S_dist)
colnames(S_dist) <- c("Iteration","Frec2008","S2008")
rownames(S_dist) <- c(1:Msim)
m <- 1
for(m in 1:Msim ){
  S_dist[m,"Iteration"] <- m
  S_dist[m,"Frec2008"] <- rpois(1, lambda=datos.frec[4,"Frecuencia"])
  S_dist[m,"S2008"] <- sum(gpdSim(model = list(xi = sev.dist[1,"xi"], mu = sev.dist[1,"mu"], beta = sev.dist[1,"beta"]), n = S_dist[m,"Frec2008"], seed = NULL))
  }

write.csv(S_dist, file = paste(path.plot,"act11302_160503_Sdist.csv", sep = ""), row.names = TRUE, col.names = TRUE)

# 		Results
sevpred.2008 <- S_dist[which(S_dist[,"S2008"]>0), ]
write.csv(sevpred.2008, file = paste(path.plot,"act11302_160503_sevpred.2008.csv", sep = ""), row.names = TRUE, col.names = TRUE)
sevhist.2008 <- hist(sevpred.2008, 50, main="Individual Severity", xlab="Simulated Data")

#
#	-- FIN: ACT11302_RuinTheory_Script.R --
