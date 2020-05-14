
# Ejemplo 

library("quantmod")

# Petrobrass series

pbr <- getSymbols("PBR", 
                  src = "yahoo", 
                  from = "2013-01-01", 
                  to = "2017-06-01", 
                  auto.assign = FALSE)
head(pbr)
class(pbr)
summary(pbr)
plot(pbr[,6])

pbr_ret <- diff(log(pbr[,6]))
pbr_ret <- pbr_ret[-1,]
plot(pbr_ret)

colnames(pbr_ret)
head(pbr_ret)

dim(pbr)
dim(pbr_ret)

# Basdaq series (USD)
spy <- getSymbols("SPY", 
                  src = "yahoo", 
                  from = "2013-01-01", 
                  to = "2017-06-01", 
                  auto.assign = FALSE)
head(spy)
class(spy)
summary(spy)
plot(spy[,6])

spy_ret <- diff(log(spy[,6]))
spy_ret <- spy_ret[-1,]
plot(spy_ret)

colnames(spy_ret)
head(spy_ret)

dim(spy)
dim(spy_ret)


# Merging returns series

port_ret <- cbind(pbr_ret,spy_ret)

head(port_ret)
dim(pbr_ret)

plot(port_ret)

summary(port_ret)


# Joint dynamics

library("mvtnorm")
library("MASS")

x.points <- seq(-.12,.12,length.out=100)
y.points <- seq(-.025,.025,length.out=100)

z <- matrix(0,nrow=100,ncol=100)

# Point estimatores
mu.est <- colMeans(port_ret)
mu.est

sigma.est <- cov(port_ret)
sigma.est

for(i in 1:100){
  for(j in 1:100){
    z[i,j] <- dmvnorm(c(x.points[i],y.points[j]),
                      mean=mu.est,
                      sigma=sigma.est)
  }
}
contour(x.points,y.points,z)


# Marginals

den.x <- dnorm(x.points, mean = mu.est[1], sd = sqrt(sigma.est[1,1]) )
den.y <- dnorm(y.points, mean = mu.est[2], sd = sqrt(sigma.est[2,2]) )


plot(y.points, den.y, 
     type="l", lty=2,
     xlab="Return value",
     ylab="Density",
     col="blue",
     main="Comparison of distributions")

lines(x.points, den.x, 
      lwd=2, col="red")

labels <- c(colnames(port_ret)[2],colnames(port_ret)[1])

legend("topright", 
       inset=.05, 
       title="Distributions",
       labels, 
       col=c("blue","red"),
       lwd=2, lty=c(1, 1))

# Portafolio

port_w <- as.matrix(c(0.65,0.35))
port_w

port_mu <- t(port_w) %*% mu.est
port_mu

port_sigma <- t(port_w) %*% sigma.est %*% port_w
port_sigma

port.points <- seq(-.1,.1,length.out=100)
den.port <- dnorm(port.points, mean = port_mu, sd = sqrt(port_sigma) )

plot(port.points, den.port, 
     type="l",
     xlab="Return value",
     ylab="Density",
     col="green",
     main="Portfolio distribution w=(.65,.36)")


