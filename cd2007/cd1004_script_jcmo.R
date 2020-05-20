
# Ejemplo 

if(!require("quantmod")){install.packages("quantmod")}
library("quantmod")

# Crude oil
crude <- getSymbols("CL=F",
                  src = "yahoo", 
                  from = "2000-03-22", 
                  to = "2016-12-31", 
                  auto.assign = FALSE)
head(crude)
class(crude)
summary(crude)
plot(crude[,6])

crude_ret <- diff(log(crude[,6]))
crude_ret <- crude_ret[-1,]
plot(crude_ret)

colnames(crude_ret)
head(crude_ret)

dim(crude)
dim(crude_ret)

# S&P500 Index

sp500 <- getSymbols("SPY", 
                  src = "yahoo", 
                  from = "2000-03-22", 
                  to = "2016-12-31", 
                  auto.assign = FALSE)
head(sp500)
class(sp500)
summary(sp500)
plot(sp500[,6])

sp500_ret <- diff(log(sp500[,6]))
sp500_ret <- sp500_ret[-1,]
plot(sp500_ret)

colnames(sp500_ret)
head(sp500_ret)

dim(sp500)
dim(sp500_ret)


# Merging returns series

port_ret <- cbind(crude_ret,sp500_ret)

head(port_ret)
dim(port_ret)

plot(port_ret)

summary(port_ret)

# Joint dynamics

if(!require("mvtnorm")){install.packages("mvtnorm")}
if(!require("MASS")){install.packages("MASS")}

library("mvtnorm")
library("MASS")

#
# Revisa los rangos de valores en summary(port_ret)
#
#> summary(port_ret)
#     Index            CL.F.Adjusted      SPY.Adjusted    
# Min.   :2000-03-23   Min.   :-0.1307   Min.   :-0.1036  
# 1st Qu.:2004-05-07   1st Qu.:-0.0129   1st Qu.:-0.0051  
# Median :2008-07-07   Median : 0.0011   Median : 0.0006  
# Mean   :2008-07-21   Mean   : 0.0004   Mean   : 0.0002  
# 3rd Qu.:2012-10-02   3rd Qu.: 0.0135   3rd Qu.: 0.0059  
# Max.   :2016-12-30   Max.   : 0.1314   Max.   : 0.1356  
#                      NA's   :1865      NA's   :945 

x.points <- seq(-.13/2,.13/2,length.out=100)
y.points <- seq(-.10/2,.13/2,length.out=100)

z <- matrix(0,nrow=100,ncol=100)

# Point estimatores
# Eliminamos los NAs (valores perdidos)
mu.est <- colMeans(port_ret, na.rm = TRUE)
mu.est

# Eliminamos los NAs (valores perdidos)
sigma.est <- cov(port_ret, use="complete.obs")
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

(sqrt(sigma.est[1,1]) > sqrt(sigma.est[2,2]))

(sqrt(sigma.est[1,1]) / sqrt(sigma.est[2,2]))


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


# CRUDE
plot(x.points, den.x, 
     type="l", lty=2,
     xlab="Return value",
     ylab="Density",
     col="blue",
     main="Petrobras")

# sp500
plot(y.points, den.y, 
     type="l", lty=2,
     xlab="Return value",
     ylab="Density",
     col="blue",
     main="Nasdaq")

# Portafolio

port_w <- as.matrix(c(.65,.35))
port_w

port_mu <- t(port_w) %*% mu.est
port_mu

port_sigma <- t(port_w) %*% sigma.est %*% port_w
port_sigma
sqrt(port_sigma)

port.points <- seq(-.1,.1,length.out=100)
den.port <- dnorm(port.points, mean = port_mu, sd = sqrt(port_sigma) )

plot(port.points, den.port, 
     type="l",
     xlab="Return value",
     ylab="Density",
     col="green",
     main="Portfolio distribution w=(.65,.35)")

# P(Ret.Port > .035) = 1 - P(Ret.Port <= .035)
1-pnorm(0.035, mean = port_mu, sd = sqrt(port_sigma))