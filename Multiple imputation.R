rm(list = ls())
require(stats4)
require(maxLik)
require(randtoolbox)
require(data.table)
require(mice)

#reading and storing data in a dataframe
dataset <- read.csv(file.choose(),header=T)
#Sample Size
N <- nrow(dataset) 
#Dependent variable (crash counts in this example); change the variable as required
DVar <- dataset$Crash 
#variable with missing data
MVar <- dataset$NLane 

#Number of imputations
M <- 100 
# index of missing values 
Mindex = which(MVar %in% NA)

#placeholders
matimpvar = matrix(NA,nrow = N, ncol = M)

# Halton Draws 
preparedraws=function()
{
  d=1
  while(d<(length(normaldraws)+1))
  {
    draws1[,normaldraws[d]]<<- qnorm(draws1[,normaldraws[d]])
    d=d+1
  }
}

Ndraws=500      # set number of draws 
dimensions=2    # define number of random parameters in the model

# generate draws (using Halton)
draws1=as.matrix(halton(Ndraws*N,dimensions))

# assign names to individual sets of draws - need one entry per dimension
colnames(draws1)=c("HRbeta1","HRbeta2")
# define whether any draws should be transformed to Normals, which is also needed for e.g. lognormals (leave empty if not)
normaldraws=c("HRbeta1","HRbeta2")

# preparing draws for estimation - this may take a while
preparedraws()

draws1 = draws1[,1:dimensions]

# prediction matrix (put 1 for those variables to be used in imputation) 
prd <- matrix(0,ncol = NCOL(dataset), nrow = NROW(dataset))
prd[5,2]=1
prd[5,4]=1

# imputation 
imp <- mice(dataset, m = M , method = "norm",predictorMatrix = prd)

# plotting imputed data vs observed parts of variable with missing data 
for (m in 1:M) {
  impdata = complete(imp,m)
  MVar2 = impdata$NLane #the variable with imputed values
  matimpvar[,m] = MVar2
}

impvar = rowMeans(matimpvar)
plot(impvar,MVar)

# plotting densities of imputed data vs observed parts of variable with missing data 
par(font=1,ps=22,family="serif",mar=c(7,7,4,2))
plot(density(MVar), col="gray",lty=1,lwd=3,
     xlab = "Number of lanes (NL)",ylab = "Density",main="",yaxt='n',ylim=c(0,0.4))
lines(density(impvar), lwd=3,col="black",lty=3)
Axis(side=2,at = seq(0,1,by=0.05))
legend(4.5,0.40, c("Density of observed variable","Density of imputed variable"),
       lty=c(1,3), lwd=c(3,3),col = c("gray","black"),bty="n",y.intersp=3)

# initial values #
init <- c(2,0.1,0.1)
#placeholders
EST = matrix(NA,nrow = NROW(init), ncol = M)
SQR = matrix(NA,nrow = NROW(init), ncol = M)
SQ = matrix(NA,nrow = NROW(init), ncol = M)

# multiple imputation
for (m in 1:M) {

# Data Preparation
impdata = complete(imp,m)
dataF =  as.matrix(data.frame(1,log(impdata$Length),impdata$NLane))

LL <- function(params){## Log Likelihood Function
    
  disp <- params[1] # Dispersion Parameter
  Fbeta <- params[2:3] # Fixed parameters in Mu Function

  offset = rep.int(dataF%*%as.matrix(Fbeta,ncol=1),Ndraws)
  mu <- exp(offset)
    
  loglik <-  sum(dnbinom(DVar,size=disp,mu=mu,log = T))
    
  return(loglik)
  }

  fit1 <- maxLik(LL,start=init,method="BFGS")
  
  EST[,m] <- fit1$estimate
  SQR[,m] <- diag(solve(-fit1$hessian))
  Lik[m] <- fit1$maximum
  
}

# Pooling estimates
PEST <- rowMeans(EST)
#Within variance
WVAR <- rowMeans(SQR)
#Between variance
for (m in 1:M) {SQ[,m] <- (EST[,m] - PEST)^2} 
BVAR <- rowSums(SQ)/(M-1)
#Total variance
TVAR <- WVAR + BVAR + BVAR/M
# Table of results
results = data.frame(PEstimates=PEST, WSD=sqrt(WVAR), BSD=sqrt(BVAR), TSD=sqrt(TVAR))

results


