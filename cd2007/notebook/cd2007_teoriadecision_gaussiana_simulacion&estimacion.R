if(!require("MASS")){install.packages("MASS")}
library("MASS")


mean.v <- rep(0.9, 0.2)
mean.v


Sigma <- matrix(c(100,-63,-63,200),2,2)
Sigma

det(Sigma)

svd(Sigma)

Precision <- solve(Sigma)
Precision

data.sim <- mvrnorm(n=3000, rep(0, 2), Sigma) 
data.sim <- as.matrix(data.sim)
data.sim <- as.data.frame(data.sim)
colnames(data.sim) <- c("Variable 1", "Variable 2")
head(data.sim)

portafolio <-  .7 * data.sim[,1] + .3 * data.sim[,2]
ts.plot(portafolio)
hist(portafolio)
media.port <- mean(portafolio)
media.port

varianza.port <- var(portafolio)
varianza.port

precision.prot <- solve(varianza.port)
precision.prot

# Propiedad 1
cport <- c(.7,.3)
cport

retorno.mean.port <- cport %*% Mean
retorno.mean.port

retorno.var.port <- t(cport) %*% Sigma %*% cport
retorno.var.port

retorno.prec.port <- solve(retorno.var.port)
retorno.prec.port

plot(data.sim)
plot(data.sim[,1])
plot(data.sim[,2])

hist(data.sim[,1])
hist(data.sim[,2])

ts.plot(data.sim[,1])
ts.plot(data.sim[,2])

# Point estimates

var.est <- var(data.sim)
var.est

lamda.est <- solve(var.est)
lamda.est

mu.est <- colMeans(data.sim)
mu.est