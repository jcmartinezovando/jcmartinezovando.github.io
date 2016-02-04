#
#	EST-46114:	Inferencia Bayesiana en Alta Dimension (Maestria en Ciencia de Datos)
#	Autor: 		Juan~Carlos Martinez Ovando
#	Email:		juan.martinez.ovando@itam.mx
#				jc.martinez.ovando@gmail.com
#	
#	Seleccion Estocastica de Variables
#	

rm(list = ls())
rutawork = 'C:/JCMO.Academia/@Cursos/2016-I_Inferencia Bayesiana en Alta Dimension/_data&code/SSVS/'
#	rutawork <- '/media/jcmo/ROCA\ ADATA/JCMO.Academia/@Cursos/2016-I_Inferencia\ Bayesiana\ en\ Alta Dimension/_data&code/SSVS/'

#	----------------------------------------------------
#		Libraries
#	----------------------------------------------------
install.packages('mvtnorm')
library(mvtnorm)

#	----------------------------------------------------
#		Datos
#	----------------------------------------------------
dde <- read.csv(paste(rutawork,'/EST46114_SSVS_Data.csv',sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)

# 	Keys
# x = dde dose
# Y = 0/1 indicator of preterm birth
# z1-z5 = possible confounders

# 	Normalization 
dde$x<- (dde$x - mean(dde$x))/sqrt(var(dde$x)) 

#	--------------------------------------------			
#			Frecuentist analysis
#	--------------------------------------------			
X <- cbind(1,dde$x,dde$z1,dde$z2,dde$z3,dde$z4,dde$z5)
Y <- dde$y
dde_mle<- glm(Y ~ -1+X, family=binomial("probit"))

# Summary table 
#Coefficients:
#            		Estimate	Std. Error 		z value		Pr(>|z|)
#(Intercept)	-1.08068	0.04355 		-24.816  	< 2e-16 ***
#dde			 0.17536	0.02909   		  6.028 		1.67e-09 ***
#z1				-0.12817	0.03528		 -3.633 	0.000280 ***
#z2				 0.11097	0.03366		  3.297 		0.000978 ***
#z3				-0.01705	0.03405		 -0.501 	0.616659
#z4				-0.08216	0.03576		 -2.298 	0.021571 *
#z5				 0.05462 	0.06473		  0.844 		0.398721

# Wald test is highly significant 
beta_mle<- dde_mle$coef # maximum likelihood estimates

#	--------------------------------------------			
#			Bayesian analysis
#	--------------------------------------------			

# 		Prior specification 
beta0 <- rep(0,7)
Pbeta0 <- 0.25*diag(7)

#		Output in beta.out
beta <- rep(0,7)	# starting value of chain
n <- nrow(dde)		# number of subjects
z <- rep(0,n)		# initial values of underlying variables
G <- 1000			# number of MCMC iterations

# 		Gibbs sampler
gtt <- 1
for(gt in 1:G){
  eta<- X%*%beta 	# linear predictor
  
  #	Sample underlying normal variables from truncated normal 
  #	from the full conditional posterior distributions
  z[Y==0] <- qnorm(runif(sum(1-Y), 0,pnorm(0,eta[Y==0],1)), eta[Y==0],1)
  z[Y==1] <- qnorm(runif(sum(Y), pnorm(0,eta[Y==1],1),1), eta[Y==1],1)

  #	Sample betas from normal full conditional posterior distribution
  Vbeta <- solve(Pbeta0 + t(X)%*%X)
  Ebeta <- Vbeta%*%(Pbeta0%*%beta0 + t(X)%*%z)
  beta <- c(rmvnorm(1,Ebeta,Vbeta))

  #	Output results
  write(t(beta),file=paste(rutawork,'/EST46114_SSVS_beta.out',sep = ""), ncol=7, append=T)
  print(c(gt,round(beta*100)/100))
}

#	Trace-plots for beta[1] and beta[2]
pdf(paste(rutawork,'/EST46114_SSVS_beta_traces.pdf',sep = ""),width=7,height=5)
par(mfrow=c(2,1))
beta_out<- matrix(scan(paste(rutawork,'/EST46114_SSVS_beta.out',sep = "")), ncol=7, byrow=T)
plot(beta_out[,1],type="l",xlab="iteration",
		ylab="intercept (beta_1)")
plot(beta_out[,2],type="l",xlab="iteration",ylab="slope (beta_2)")
dev.off()

#	Plot of marginal posterior density of slope 
pdf(paste(rutawork,'/EST46114_SSVS_beta_marginal_slope.pdf',sep = ""),width=7,height=5)
slp <- beta_out[101:1000,2]
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
plot(density(slp),type="l",xlab="DDE slope (beta_1)",ylab="Posterior Density",
     cex=1.2)
abline(v=mean(slp))
abline(v=0.17536139,lwd=2.5) # MLE
abline(v=0.17536139 + 1.96*0.02909*c(-1,1),lwd=2.5,lty=2)
abline(v=quantile(slp,probs=c(0.025,0.975)),lty=2)
dev.off()

#	Plot estimated dose response curve for preterm birth
pdf(paste(rutawork,'/EST46114_SSVS_beta_doseresponse.pdf',sep = ""),width=7,height=5)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
xg<- seq(-2,2,length=100)  # grid of dde values
beta<- beta_out[101:1000,] # discard burn-in
post<- matrix(0,100,4)
for(i in 1:100){
  post[i,1]<- mean(pnorm(beta[,1] + xg[i]*beta[,2]))
  post[i,2:3]<- quantile(pnorm(beta[,1] + xg[i]*beta[,2]),probs=c(0.025,0.975))
  post[i,4]<- pnorm(-1.08068 + xg[i]*0.17536139)
}
xtrue<- xg*sqrt(var(dde$x)) + mean(dde$x) # back transform
plot(xtrue,post[,1],xlab="Serum DDE (mg/L)",ylab="Pr(Preterm birth)",cex=1.2,
     ylim=c(0,max(post)), type="l")
lines(xtrue,post[,2],lty=2)
lines(xtrue,post[,3],lty=2)
lines(xtrue,post[,4],lwd=2.5)
dev.off()

#	Posterior summaries of regression coefficients
table1<- matrix(0,7,5)
for(i in 1:7){
  table1[i,]<- c(mean(beta[,i]),median(beta[,i]),sqrt(var(beta[,i])),
                 quantile(beta[,i],probs=c(0.025,0.975)))
}
table1<- round(table1*100)/100
write.csv(t(table1),file=paste(rutawork,'/EST46114_SSVS_beta_summary.csv',sep = ""))

#	---------------------------------------------------------------------			
#		Stochastic search variable selection via the Gibbs sampler
#	---------------------------------------------------------------------			

ddemle<- glm(Y ~ -1+X, family=binomial("probit"))
betamle<- ddemle$coef 		#	MLE

#	Prior specification
p <- ncol(X)        # number of predictors
p0 <- rep(0.5,p)    # prior probability of excluding a predictor
b0 <- rep(0,p)      # prior mean for normal component if predictor included
s0 <- rep(2,p)      # standard deviation for normal component
	
#	Stochastic search variable selection via data-augmentation
beta <- betamle
n <- nrow(dde)   	# number of subjects
z <- rep(0,n)    	# initial values of underlying variables
G <- 5000        	# number of MCMC iterations

# -------------------------------------
#	Stochastic search Gibbs sampler
# -------------------------------------
gt <- 1
for(gt in 1:G){
  #		Data augmentation step
  eta<- X%*%beta 	# linear predictor
  
  # Sample underlying normal variables from truncated normal 
  # from the full conditional posterior distributions
  z[Y==0]<- qnorm(runif(sum(1-Y),0,pnorm(0,eta[Y==0],1)),eta[Y==0],1)
  z[Y==1]<- qnorm(runif(sum(Y),pnorm(0,eta[Y==1],1),1),eta[Y==1],1)

  #	Update coefficients one at a time (unlike in the previous algorithm)
  j <-1
  for(j in 1:p){
    #	Conditional posterior variance of beta_j under normal prior
    Vj<- 1/(s0[j]^{-2} + sum(X[,j]^2))
    Ej<- Vj*sum(X[,j]*(z-X[,-j]%*%beta[-j]))
    
	#	Conditional probability of including jth predictor
    phat<- 1/(1 + p0[j]/(1-p0[j]) * dnorm(0,Ej,sqrt(Vj))/dnorm(0,b0[j],s0[j]) )          
    m<- rbinom(1,1,phat)
    beta[j]<- m*rnorm(1,Ej,sqrt(Vj))
  }
  
  #	Output results
  write(t(beta),file=paste(rutawork,'/EST46114_SSVS_beta_ss.out',sep = ""), ncol=7, append=T)
  print(c(gt,round(beta*100)/100))
}

#	Trace-plots for beta[1] and beta[2]
beta<- matrix(scan(paste(rutawork,'/EST46114_SSVS_beta_ss.out',sep = "")), ncol=7, byrow=T)
# 2^p possible models -> all possible combinations of excluding each 
# of the p predictors, including the intercept 
# In this case, 2^p = 2^7 = 128

pdf(paste(rutawork,'/EST46114_SSVS_beta_traces_ss.pdf',sep = ""),width=7,height=5)
# Plot Gibbs iterations 
par(mfrow=c(3,2))
ylb=c("intercept","dde","z1","z2","z3","z4","z5")
for(j in 2:7){
  print(j)
  plot(beta[,j],xlab="iteration",ylab=ylb[j])
}
dev.off()

#	Plot estimated dose response curve for preterm birth
pdf(paste(rutawork,'/EST46114_SSVS_beta_doseresponse_ss.pdf',sep = ""),width=7,height=5)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
xg<- seq(-2,2,length=100)  # grid of dde values
beta<- beta[1001:nrow(beta),] # discard burn-in
post<- matrix(0,100,4)
for(i in 1:100){
  post[i,1]<- mean(pnorm(beta[,1] + xg[i]*beta[,2]))
  post[i,2:3]<- quantile(pnorm(beta[,1] + xg[i]*beta[,2]),probs=c(0.025,0.975))
  post[i,4]<- pnorm(-1.08068 + xg[i]*0.17536139)
}
xtrue<- xg*sqrt(var(dde$x)) + mean(dde$x) # back transform
plot(xtrue,post[,1],xlab="Serum DDE (mg/L)",ylab="Pr(Preterm birth)",cex=1.2,
     ylim=c(0,max(post)), type="l")
lines(xtrue,post[,2],lty=2)
lines(xtrue,post[,3],lty=2)
lines(xtrue,post[,4],lwd=2.5)
dev.off()

#	Posterior summaries of regression coefficients
#	Columns include:
#		- posterior mean
#		- median
#		- standard deviation
#		- 95% credible interval
#		- Pr(beta_j=0|data)
table1<- matrix(0,7,6)
for(i in 1:7){
  table1[i,]<- c(mean(beta[,i]),median(beta[,i]),sqrt(var(beta[,i])),
                 quantile(beta[,i],probs=c(0.025,0.975)),
		 length(beta[beta[,i]==0,i])/nrow(beta))
}
table1<- round(table1*100)/100
write(t(table1),file=paste(rutawork,'/EST46114_SSVS_beta_summary_ss.csv',sep = ""),ncol=6)

#	Calculate posterior probabilities for the best models among those visited
Mout <- beta 
Mout <- matrix(as.numeric(I(Mout==0)),nrow(Mout),ncol(Mout)) 
Mindex <- Mout[1,] 		# different models sampled starting with first
						# returns 1 if m1 and m2 are the same model
ind<- function(m1,m2){
  as.numeric(all(I(m1==m2)))
}
Im<- apply(Mout,1,ind,m2=Mout[1,]) # indicators of samples from 1st model
Nm<- sum(Im)                       # number of samples from 1st model
Mout<- Mout[Im==0,]
repeat{
  if(length(Mout)==7) Mout<- matrix(Mout,1,7)
  Mindex<- rbind(Mindex,Mout[1,])
  Im<- apply(Mout,1,ind,m2=Mout[1,])
  Nm<- c(Nm,sum(Im))
  print(sum(Nm))
  if(sum(Nm)==nrow(beta)){ 
    break
  } else Mout<- Mout[Im==0,]
}
#	Sort models visited by decreasing posterior probability
Pm <- Nm/sum(Nm)
ord <- order(Pm)
Pm <- Pm[rev(ord)]
Mindex <- Mindex[rev(ord),]
table2 <- cbind(Pm,Mindex)
write(t(table2),file=paste(rutawork,'/EST46114_SSVS_beta_pmodel_ss.csv',sep = ""),ncol=8)

#
#	-- END: EST46114_SSVS_Script.r --