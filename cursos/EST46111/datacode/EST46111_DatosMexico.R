#
#	EST-46111:	Fundamentos de Estadística (Maestría en Ciencia de Datos)
#	Autor: 			Juan~Carlos Martínez Ovando
#	Email:			juan.martinez.ovando@itam.mx
#						jc.martinez.ovando@gmail.com
#	
#	Modelos Lineale de Regresión (Enfoque Bayesiano)
#	

rm(list = ls())
rutawork = 'C:/JCMO.Academia/@Cursos/2015-II_Fundamentos Estadistica/_data&code/Mexico'
 
#	----------------------------------------------------
#		Código
#	----------------------------------------------------
source(paste(rutawork,"/baylinreg.R",sep = ""))
library("MASS")

#	----------------------------------------------------
#		Datos
#	----------------------------------------------------
#	Loading the data...
rdata <- read.csv(paste(rutawork,"/EST46111_DatosMexico.csv",sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)

rdata <- as.data.frame(rdata)

plot(rdata$inpc_general, rdata$salario)
plot(rdata$mexindpro, rdata$salario)

Y <- as.matrix(rdata$salario)
X <- as.matrix(cbind(matrix(1,dim(Y)), rdata$inpc_general,rdata$mexindpro))

colnames(Y) <- c("salarios")
colnames(X) <- c("cte","inpc","prodind")

#	Initialization
T <- dim(Y)[1]
p <- dim(X)[2]

M <- 1000 			#	Número de simulaciones

#	----------------------------------------------------
#		Regresión lineal simple
#	----------------------------------------------------
m_0 <-  solve(t(X)%*%X)%*%t(X)%*%Y
S_0 <- diag(1,p,p)
a_0 <- 1
b_0 <- 1
salarios.linreg <- baylinreg(Y,X,m_0,S_0,a_0,b_0)

lambda_sim <- rgamma(M, shape=salarios.linreg[[3]], scale = 1/salarios.linreg[[4]])
m_sim <- matrix(NaN,M,p)
t <- 1
for(t in 1:M){
	m_sim[t,] <- mvrnorm(n = 1, mu=salarios.linreg[[1]], Sigma=solve(lambda_sim[t]*salarios.linreg[[2]]))
  }
colnames(m_sim) <- colnames(X)  

#	----------------------------------------------------
#		Regresión lineal con efectos fijos
#	----------------------------------------------------
J <- max(unique(rdata[,"clave"])) - 1
rdataf <- rdata
j <-1 
for(j in 1:J){
	f <- 0*rdata[,1]
	f[which(rdata[,"clave"]==j)] <- 1
	rdataf <- cbind(rdataf,f)
  }
Xf <- as.matrix(cbind(matrix(1,dim(Y)), rdataf$inpc_general,rdataf$mexindpro,rdataf[,c((6+1):(6+J))]))
  
mf_0 <-  solve(t(Xf)%*%Xf)%*%t(Xf)%*%Y
pf <- ncol(Xf)
Sf_0 <- diag(1,pf,pf)
af_0 <- 1
bf_0 <- 1
salarios.linreg.fixed <- baylinreg(Y,Xf,mf_0,Sf_0,af_0,bf_0)

lambdaf_sim <- rgamma(M, shape=salarios.linreg.fixed[[3]], scale = 1/salarios.linreg.fixed[[4]])
mf_sim <- matrix(NaN,M,pf)
t <- 1
for(t in 1:M){
	mf_sim[t,] <- mvrnorm(n = 1, mu=salarios.linreg.fixed[[1]], Sigma=solve(lambdaf_sim[t]*salarios.linreg.fixed[[2]]))
  }
colnames(mf_sim) <- colnames(Xf)  

hist(m_sim[,2])
hist(mf_sim[,2])

#	----------------------------------------------------
#		Salidas
#	----------------------------------------------------
save( salarios.linreg,
		 m_sim, lambda_sim,
		 salarios.linreg.fixed,
		 mf_sim, lambdaf_sim,
	     file = paste(rutawork,"/salarios_output.RData",sep = "")
	)
	
#
#		FIN de "EST46111_DatosMexico.R"