#
#	This program fits a hierarchical dynamic bayesian model averaging. 
#	The output variable is the States' economic growth for Mexico 
#	for the 32 Mexican States. Output is measured using the ITAEE data
#	in forward annual variatios on a quarterly basis.
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez.ovando@itam.mx
#				JC.Martinez.Ovando@gmail.com
#
#	Reference:	Martínez-Ovando, J.~C. (2014) "The Influence of Crime on the Economic Growth 
#					Dynamics in Mexico in the Short-Term," Banco de México, Mimeo.
#
#	Notes:		-	Data must be processed before loading it in R 
#					(e.g. growth rates, truncation periods, imputations, etc.)
#				-	Data spans from 2006Q4 to 2011Q3 (i.e. 20 quarterly observations)
#				
#		

rm(list = ls())
rutawork = 'C:/JCMO.Trabajo/@Mis.Cursos/2016-I_Inferencia Bayesiana en Alta Dimension/_data&code/dlmbma'
rutacode = 'C:/JCMO.Trabajo/@Mis.Cursos/2016-I_Inferencia Bayesiana en Alta Dimension/_data&code/dlmbma/code'
 
#	----------------------------------------------------
#		Loading code
#	----------------------------------------------------
source(paste(rutacode,"/dlm.bma.nest.R",sep = ""))
source(paste(rutacode,"/forward.filter.R",sep = ""))
source(paste(rutacode,"/dlm.bma.post.filter.R",sep = ""))
 
#	----------------------------------------------------
#		Loading data
#	----------------------------------------------------
#	Loading the data...
rdata <- read.csv(paste(rutawork,"/input/EST46114_rdata_160519.csv",sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)

#	Loading matrix of seasonal covariates
Smat <- read.csv(paste(rutawork,"/input/EST46114_Smat_seasonal_160519.csv",sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE, row.names = 1)
Smat <- as.matrix(Smat)

Qini <- 20
Qnum <- 20
periodosana <- 	c("2006Q4",
				  "2007Q1","2007Q2","2007Q3","2007Q4",
				  "2008Q1","2008Q2","2008Q3","2008Q4",
				  "2009Q1","2009Q2","2009Q3","2009Q4",
				  "2010Q1","2010Q2","2010Q3","2010Q4",
				  "2011Q1","2011Q2","2011Q3")
  
#	----------------------------------------------------
#		Repository for the output of 32 States
#	----------------------------------------------------

#	Object lists
dma.Xs <- list()

dlmbma.y.filt <- list()
dlmbma.Q.filt <- list()
dlmbma.y.fits <- list()
dlmbma.yQ.fits <- list()
dlmbma.pi.filt <- list()
dlmbma.data.j <- list()
dlmbma.mt.ma <- list()
dlmbma.mt.ma.one <- list()
dlmbma.mt.ma.two <- list()

#	Correlation of adjusted values (MLS estimation, no inference)
rho <- matrix(NA,32,2)
colnames(rho) <- c("Adj.Corr.Y.one","Adj.Corr.Y.two")
rownames(rho) <- c(1:32)

#	----------------------------------------------------
#		Running the model	-	State level
#	----------------------------------------------------
J <- 32 					#	Number of Mexican States
ent <- 1
for(ent in 1:J){
		#'''''''''''''''''''''''''''''''''''''''''''''''		
		#		Setting-up the data set
		#'''''''''''''''''''''''''''''''''''''''''''''''

		#		Indexes
		indexes <- which(rdata[,"clave"] == ent)
		
		#		Annual growth rate (a period ahead)
		Y <- rdata[indexes,c("periodo","itaee_c")]
		Ydat <- as.matrix(Y[Qini:(Qini+Qnum-1),"itaee_c"])
		colnames(Ydat) <- c("itaee_c")
		rownames(Ydat) <- periodosana
		
		#		Covariates X's (initial conditions)
		X <- rdata[indexes,]
		Xmat <- X[Qini:(Qini+Qnum-1), c('itaee_log','itaee_c1','itaee_c2','inpc_general_c1','desocup','ied_mxp_log','usindpro_log','pibp12mes_e','gtopgj_pc')]
		rownames(Xmat) <- periodosana
	
		#		Design matrix	(Note: the function "dyn.dma.R" automatically introduces a constant variable, before S´s, X´s and C´s variables)
		Xdat <- cbind(Smat,Xmat)

		Data_j <- cbind(Ydat,Xdat)
		
		Xdat <- as.matrix(Xdat)
		dma.Xs[[ent]] <- Xdat
		Best <- solve(t(Xdat)%*%Xdat)%*%t(Xdat)%*%Ydat
		Yhat.one <- Xdat %*%solve(t(Xdat)%*%Xdat)%*%t(Xdat)%*%Ydat								#	with CRIME variables	
		Yhat.two <- Xdat[,1:13] %*%solve(t(Xdat[,1:13])%*%Xdat[,1:13])%*%t(Xdat[,1:13])%*%Ydat	#	without CRIME variables
		
		rho[ent,1] <- cor(Ydat,Yhat.one) 
		rho[ent,2] <- cor(Ydat,Yhat.two) 

		#	--------------------------------------------
		#		14-09-15	-	Re defining the selection matrix	
		#	--------------------------------------------
		#		Defining which covariables to switch from
		Xwhich <- read.csv(paste(rutawork,"/input/EST46114_Xwhich_160519.csv",sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
		modelthreshold <- 21

		#	Executing Dynamic Bayesian Model Averaging
		x <- Xdat[,-1]
		y <- Ydat
		models.which <- Xwhich 
		lambda <- 0.999
		gammaa <- 0.999
		eps <- 0.001/nrow(models.which)

		#'''''''''''''''''''''''''''''''''
		#	Function:	dlm.bma.nest
		#'''''''''''''''''''''''''''''''''
		dlmbma_ent <- dlm.bma.nest(x, y, models.which, lambda, gammaa, eps, modelthreshold)
		
		#	Saving the output
		#	Filtered
		y.fits <- cbind(Ydat,Yhat.one,Yhat.two,dlmbma_ent$y.filt.ma,dlmbma_ent$y.filt.ma.one,dlmbma_ent$y.filt.ma.two)
		colnames(y.fits) <- c("itaee_c","itaee_yhat.one","itaee_yhat.two","itaee_yfilt_ma","itaee_yfilt_ma.one","itaee_yfilt_ma.two")

		#	Predictive dispersion
		yQ.fits <- cbind(Ydat,dlmbma_ent$y.filt.ma,sqrt(dlmbma_ent$Q.filt.ma),dlmbma_ent$y.filt.ma.one,sqrt(dlmbma_ent$Q.filt.ma.one),dlmbma_ent$y.filt.ma.two,sqrt(dlmbma_ent$Q.filt.ma.two))
		colnames(yQ.fits) <- c("itaee_c","itaee_yfilt_ma","itaee_Qfilt_ma","itaee_yfilt_ma.one","itaee_Qfilt_ma.one","itaee_yfilt_ma.two","itaee_Qfilt_ma.two")
		
		dlmbma.y.filt[[ent]] <- dlmbma_ent$y.filt
		dlmbma.Q.filt[[ent]] <- dlmbma_ent$Q.filt
		dlmbma.y.fits[[ent]] <- y.fits
		dlmbma.yQ.fits[[ent]] <- yQ.fits
		dlmbma.pi.filt[[ent]] <- dlmbma_ent$pi.filt
		dlmbma.data.j[[ent]] <- Data_j
		dlmbma.mt.ma[[ent]] <- dlmbma_ent$m.filt.ma
		dlmbma.mt.ma.one[[ent]] <- dlmbma_ent$m.filt.ma.one
		dlmbma.mt.ma.two[[ent]] <- dlmbma_ent$m.filt.ma.two
		
	#	 END of dlm.bma for the 32 States
	}
	
#''''''''''''''''''''''''''''''''''''''''''''''''
#	Exporting CSV outputs
#''''''''''''''''''''''''''''''''''''''''''''''''
#	Exporting the mt's for general BMA and scenarios ONE and TWO
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.mt.ma[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output/mt_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	matrix.aux <- dlmbma.mt.ma.one[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output/mt_",ent,"_one.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	matrix.aux <- dlmbma.mt.ma.two[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output/mt_",ent,"_two.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
#	Exporting Data_j
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.data.j[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output/data_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
#	Exporting posterior filtering paths
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.y.filt[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output/filter_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	matrix.aux <- dlmbma.Q.filt[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output/filter_",ent,"_Q.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
#	Exporting posterior filtering paths
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.y.fits[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output/fit_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	matrix.aux <- dlmbma.yQ.fits[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output/fit_",ent,"_Q.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
#	Exporting posterior probabilities for the alternative models
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.pi.filt[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output/pmp_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
	
#	Saving outputs
save( dlmbma.y.filt, dlmbma.Q.filt,
	  dlmbma.y.fits, dlmbma.yQ.fits,
	  dlmbma.pi.filt,
	  dlmbma.data.j,
	  dlmbma.mt.ma,
	  dlmbma.mt.ma.one,
	  dlmbma.mt.ma.two,
	  file = paste(rutawork,"/output/EST46114_dlmbma_160519.RData",sep = "")
	)
	
#
#		END of "EST46114_dlmbma_160519.R"