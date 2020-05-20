
# series para elegir: DOS de ellas 

if(!require("quantmod")){install.packages("quantmode")}
library("quantmod")

# Gold June 20
gold <- getSymbols("GC=F",
                  src = "yahoo", 
                  from = "2000-03-23", 
                  to = "2016-12-30", 
                  auto.assign = FALSE,)
plot(gold[,6])
summary(gold)

# Crude oil
crude <- getSymbols("CL=F",
                  src = "yahoo", 
                  from = "2000-03-23", 
                  to = "2016-12-30", 
                  auto.assign = FALSE)

plot(crude[,6])
summary(crude)

# S&P500 Index
sp500 <- getSymbols("SPY", 
                  src = "yahoo", 
                  from = "2000-03-22", 
                  to = "2016-12-31", 
                  auto.assign = FALSE)
plot(sp500[,6])
summary(sp500)

# GBP/USD (GBPUSD=X)
gbpusd <- getSymbols("GBPUSD=X", 
                  src = "yahoo", 
                  from = "2000-03-22", 
                  to = "2016-12-31", 
                  auto.assign = FALSE)
plot(gbpusd[,6])
summary(gbpusd)

# FTSE 100 (^FTSE)
ftse100 <- getSymbols("^FTSE", 
                  src = "yahoo", 
                  from = "2000-03-22", 
                  to = "2016-12-31", 
                  auto.assign = FALSE)
plot(ftse100[,6])
summary(ftse100)
