##############################
## BOOTSTRAP (NONPARAMETRIC BLOCK BOOTSTRAP WITH cenred and rescaled)
## 2. Internal bootstrap: Internal Centered and Rescaled Residuals + External treatments

rm(list=ls())

library (nlme)

generic.group.boot.Internal.ExtSecondCalbraion.plus.Adjustment<- function(y,X,alpha,beta,sigma.u,sigma.e,Var.sigma.u,Var.sigma.e,Cov.simga.ue, group,k=1000,verbose=F,stop=F) {
	N <- length(group)   				# group is region or domain, N is sample size
	group.labels <- unique(group)			#unique  domain code
	H <- length(group.labels)			#total number of domains
	P <- ncol(X)
	y.fit <- alpha + X%*%as.matrix(beta,ncol=1) 	 #fitted y values
	y.resid <- y-y.fit 					 #Raw residuals
	y.resid.u1 <- rep(0,H)					#Domain level residuals
	for(h in 1:H) y.resid.u1[h] <- mean(y.resid[group==group.labels[h]])
	y.resid.e1 <- rep(0,N)
	for(h in 1:H) y.resid.e1[group==group.labels[h]] <- y.resid[group==group.labels[h]]-y.resid.u1[h]
	a.rescal=sqrt((sigma.u^2)/((t(y.resid.u1)%*%y.resid.u1)/H))
	y.resid.u=a.rescal*y.resid.u1
	b.rescal=sqrt((sigma.e^2)/((t(y.resid.e1)%*%y.resid.e1)/N))
	y.resid.e=b.rescal*y.resid.e1
	lambda <- (sigma.u/sigma.e)^2
# This matrix kept random.effects for all domains
# This matrix kept the estimated parametrs: p+1 beta ests, sigma u, sigm e and their ratio
	boot.mat.random.effects <- matrix(0,nrow=k,ncol=H)
	boot.mat.model.parameters <- matrix(0,nrow=k,ncol=P+4)


for(b in 1:k) {
## bootstrap simulations start
		y.b <- rep(0,N)
		donor.group.labels.b <- sample(x=group.labels,size=H,replace=T)
		y.resid.u.b <- sample(y.resid.u,size=H,replace=T)
		corr.fact=mean(y.resid.u)           ## Defined mean of area residual means
		y.resid.u.b=y.resid.u.b-corr.fact  #Centered Domain level residuals

for(h in 1:H) {
			target.units <- (1:N)[group==group.labels[h]]
			donor.units <- (1:N)[group==donor.group.labels.b[h]]
			if(length(donor.units)>1) donating.units <- sample(x=donor.units,size=length(target.units),replace=T)
			else donating.units. <- rep(donor.units,length(target.units))
			y.b[target.units] <- y.fit[target.units]+y.resid.u.b[h]+y.resid.e[donating.units]
		}


y.model.b <- lme(y.b~X,random = ~1 | group)
		var.comp.b <- as.real(VarCorr(y.model.b))[1:2]
		lambda.b <- var.comp.b[1]/var.comp.b[2]
		boot.mat.model.parameters[b,] <- c(y.model.b$coef$fixed,var.comp.b,lambda.b)		# This matrix kept the estimated parametrs: p+1 beta ests, sigma u, sigm e and their ratio
		boot.mat.random.effects[b,] <- c(y.model.b$coef$random$group) 			 	# This matrix kept random.effects for all domains
		if(verbose) print(c("bootstrap simulation ",b))
## bootstrap simulations end k times
	}


## Second order correction : Cholesky decomposition
Var.Sigma<-matrix(c(Var.sigma.u,Cov.simga.ue,Cov.simga.ue,Var.sigma.e),2,2)
chol.de.Var.Sigma=t(chol(Var.Sigma))
#chol.de.Var.Sigma%*%t(chol.de.Var.Sigma);Var.Sigma
boot.mat.model.adj.parameters 	<- boot.mat.model.parameters
boot.mat.model.adj.parameters1 	<- boot.mat.model.parameters
boot.mat.model.adj.parameters2	<- boot.mat.model.adj.parameters1[,(P+2):(P+3)]


## centred to zero mean
boot.mat.model.adj.parameters3	<-cbind(boot.mat.model.adj.parameters2[,1]-(apply(boot.mat.model.adj.parameters2,2,mean))[1], boot.mat.model.adj.parameters2[,2]-(apply(boot.mat.model.adj.parameters2,2,mean))[2])
Var.Sigma.boot<-matrix(c(var(boot.mat.model.adj.parameters3[,1]),cov(boot.mat.model.adj.parameters3[,1],boot.mat.model.adj.parameters3[,2]),
                         cov(boot.mat.model.adj.parameters3[,1],boot.mat.model.adj.parameters3[,2]),var(boot.mat.model.adj.parameters3[,2])),2,2)

chol.de.Var.Sigma.boot=t(chol(Var.Sigma.boot))
#chol.de.Var.Sigma.boot%*%t(chol.de.Var.Sigma.boot);Var.Sigma.boot
A=t(chol.de.Var.Sigma%*%solve(chol.de.Var.Sigma.boot))

boot.mat.model.adj.parameters4<-boot.mat.model.adj.parameters3%*%A

boot.mat.model.adj.parameters[,(P+2):(P+3)]<- cbind(boot.mat.model.adj.parameters4[,1]+(apply(boot.mat.model.adj.parameters2,2,mean))[1],boot.mat.model.adj.parameters4[,2]+(apply(boot.mat.model.adj.parameters2,2,mean))[2])

if(stop) browser()
	list(
	bootstrap.par=boot.mat.model.parameters,
	bootstrap.adj.par=boot.mat.model.adj.parameters,
	bootstrap.random=boot.mat.random.effects)

}

### Start of simulation
library (nlme)
#set.seed(100)
NoSim=10
NoBootSim=100
areaNo=20
areasize=50

test.x <- rnorm(areaNo*areasize,1,2)            # x matrix to size n=1000
test.g <- sample(rep(1:areaNo,areasize))     # Number of domains D=100
group.boot.sim1.2 <- matrix(0,nrow=NoSim,ncol=25)
mat.Var.Sigma <- matrix(0,nrow=NoSim,ncol=9)

## Boostrap mean of sigma u and e.
Bootstap.estimate.sim<- matrix(0,nrow=NoSim,ncol=4)

## for W50x50.Spatial.New
W<- as.matrix(read.table("C:/Users/Payam/Documents/W.txt"))

Rho=0
True.beta0= 1
True.beta1= 2
True.Sigma2u=1
True.Sigma2e=9
TruesigmaRatio=True.Sigma2u/True.Sigma2e


## 0.04/(0.04+True.Sigma2e)  0.1001821
##True.Sigma2e  0.3592728


sim.start <- 1
sim.finish <- NoSim

for(i in sim.start:sim.finish) {
cat(date(),"Monte Carlo iteration",i,"\n")
	test.u <- rnorm(areaNo,0,sqrt(True.Sigma2u))           				# area effects
	test.e <- rnorm(areasize*areaNo,0,sqrt(True.Sigma2e))				# indivd effects
	u.y <- solve(diag(areaNo)-Rho*W)%*%test.u
	test.z <- kronecker(u.y,rep(1,areasize))         			# replicated area effects by d times


	#mean(test.e);mean(e.y)
	#var(test.e);var(e.y)
	test.y <- 1+2*test.x+test.z+test.e 		 # y data generated

#solve((diag(areasize)-Rho*W)%*%(diag(areasize)-Rho*t(W)))


	test.fit <- lme(test.y~test.x,random = ~1 | test.g)	# model fit for sample data generated
	var.comp <- as.real(VarCorr(test.fit))[1:2]
	lambda <- var.comp[1]/var.comp[2]
	estimates <- c(test.fit$coeff$fixed[1:2],var.comp,lambda)



	#diag(test.fit$varFix)	## sd of beta estimates
	#matrix(as.real(test.fit$apVar),2,2) ## Var of varaince componetes
	est.rand.eff=c(test.fit$coef$random$test.g)


ni <- rep(0,areaNo)
for(tt in 1:areaNo)ni[tt]<- sum(test.g==tt)

vi.inv <- function(i) {solve(estimates[4]*(diag(rep(1,ni[i])))+estimates[3]*(matrix(1,ni[i],ni[i]))) }

T1=0
for (ii in 1:areaNo){
T=sum(diag((vi.inv(ii)%*%matrix(1,ni[ii],ni[ii]))%*%vi.inv(ii)%*%matrix(1,ni[ii],ni[ii])))
T1=T1+T
}

T2=0
for (ii in 1:areaNo){
T=sum(diag(vi.inv(ii)%*%vi.inv(ii)))
T2=T2+T
}

T12=0
for (ii in 1:areaNo){
T=sum(diag((vi.inv(ii)%*%matrix(1,ni[ii],ni[ii]))%*%vi.inv(ii)))
T12=T12+T
}

Ju<-(1/2)*T1
Je<-(1/2)*T2
Jue<-(1/2)*T12
DD<-Ju*Je-(Jue^2)

# Var-Cov of estimates of variance components
Var.u.hat  <- (1/DD)*Je                        # Var(estsigma2u)
Var.e.hat  <- (1/DD)*Ju	                  	 # Var(estsigma2e)
Var.ue.hat <-(-1/DD)*Jue 	                   # Cov(estsigma2u, estsigma2e)

mat.Var.Sigma[i,1] <-Var.u.hat
mat.Var.Sigma[i,2] <-Var.e.hat
mat.Var.Sigma[i,3] <-Var.ue.hat


	test.boot1.2 <- generic.group.boot.Internal.ExtSecondCalbraion.plus.Adjustment(y=test.y,X=as.matrix(test.x,ncol=1),alpha=estimates[1],beta=estimates[2],sigma.u=sqrt(estimates[3]),sigma.e=sqrt(estimates[4]),Var.sigma.u=Var.u.hat,Var.sigma.e=Var.e.hat,Cov.simga.ue=Var.ue.hat,group=test.g,k=NoBootSim)

group.boot.sim1.2[i,1:5] <- estimates         ## it sample estimates of the parametrs  fixed effect,  sigma.u, sigma.e, their ratio

group.boot.sim1.2[i,6]  <-(quantile(test.boot1.2$bootstrap.par[,1],probs=0.025)<True.beta0)*(quantile(test.boot1.2$bootstrap.par[,1],probs=0.975)>True.beta0) 	 	##CR of beta0
group.boot.sim1.2[i,7] <- (quantile(test.boot1.2$bootstrap.par[,2],probs=0.025)<True.beta1)*(quantile(test.boot1.2$bootstrap.par[,2],probs=0.975)>True.beta1)  		##CR of beta1
group.boot.sim1.2[i,8] <- (quantile(test.boot1.2$bootstrap.par[,3],probs=0.025)<True.Sigma2u)*(quantile(test.boot1.2$bootstrap.par[,3],probs=0.975)>True.Sigma2u)	##CR of sigmau
group.boot.sim1.2[i,9] <- (quantile(test.boot1.2$bootstrap.par[,4],probs=0.025)<True.Sigma2e)*(quantile(test.boot1.2$bootstrap.par[,4],probs=0.975)>True.Sigma2e)	##CR of sigmae
group.boot.sim1.2[i,10] <-(quantile(test.boot1.2$bootstrap.par[,5],probs=0.025)<TruesigmaRatio)*(quantile(test.boot1.2$bootstrap.par[,5],probs=0.975)>TruesigmaRatio)#CR embda


## CR for adjusted parameters estimates
	group.boot.sim1.2[i,11] <- (quantile(test.boot1.2$bootstrap.adj.par[,1],probs=0.025)<True.beta0)*(quantile(test.boot1.2$bootstrap.adj.par[,1],probs=0.975)>True.beta0)
	group.boot.sim1.2[i,12] <- (quantile(test.boot1.2$bootstrap.adj.par[,2],probs=0.025)<True.beta1)*(quantile(test.boot1.2$bootstrap.adj.par[,2],probs=0.975)>True.beta1)
	group.boot.sim1.2[i,13] <- (quantile(test.boot1.2$bootstrap.adj.par[,3],probs=0.025)<True.Sigma2u)*(quantile(test.boot1.2$bootstrap.adj.par[,3],probs=0.975)>True.Sigma2u)
	group.boot.sim1.2[i,14] <- (quantile(test.boot1.2$bootstrap.adj.par[,4],probs=0.025)<True.Sigma2e)*(quantile(test.boot1.2$bootstrap.adj.par[,4],probs=0.975)>True.Sigma2e)
	group.boot.sim1.2[i,15] <- (quantile(test.boot1.2$bootstrap.adj.par[,5],probs=0.025)<TruesigmaRatio)*(quantile(test.boot1.2$bootstrap.adj.par[,5],probs=0.975)>TruesigmaRatio)

## width for parameters estimates
	group.boot.sim1.2[i,16]  <- quantile(test.boot1.2$bootstrap.par[,1],probs=0.975)-quantile(test.boot1.2$bootstrap.par[,1],probs=0.025) ## width  of beta0
	group.boot.sim1.2[i,17]  <- quantile(test.boot1.2$bootstrap.par[,2],probs=0.975)-quantile(test.boot1.2$bootstrap.par[,2],probs=0.025) ## width of beta1
	group.boot.sim1.2[i,18]  <- quantile(test.boot1.2$bootstrap.par[,3],probs=0.975)-quantile(test.boot1.2$bootstrap.par[,3],probs=0.025) ## width of sigmau
	group.boot.sim1.2[i,19]  <- quantile(test.boot1.2$bootstrap.par[,4],probs=0.975)-quantile(test.boot1.2$bootstrap.par[,4],probs=0.025) ## width of sigmae
	group.boot.sim1.2[i,20]  <- quantile(test.boot1.2$bootstrap.par[,5],probs=0.975)-quantile(test.boot1.2$bootstrap.par[,5],probs=0.025) ## width of lembda

## width for adjusted parameters estimates
	group.boot.sim1.2[i,21]  <- quantile(test.boot1.2$bootstrap.adj.par[,1],probs=0.975)-quantile(test.boot1.2$bootstrap.adj.par[,1],probs=0.025) ## width  of adjusted beta0
	group.boot.sim1.2[i,22]  <- quantile(test.boot1.2$bootstrap.adj.par[,2],probs=0.975)-quantile(test.boot1.2$bootstrap.adj.par[,2],probs=0.025) ## width of  adjusted beta1
	group.boot.sim1.2[i,23]  <- quantile(test.boot1.2$bootstrap.adj.par[,3],probs=0.975)-quantile(test.boot1.2$bootstrap.adj.par[,3],probs=0.025) ## width of adjusted sigmau
	group.boot.sim1.2[i,24]  <- quantile(test.boot1.2$bootstrap.adj.par[,4],probs=0.975)-quantile(test.boot1.2$bootstrap.adj.par[,4],probs=0.025) ## width of adjusted sigmae
	group.boot.sim1.2[i,25]  <- quantile(test.boot1.2$bootstrap.adj.par[,5],probs=0.975)-quantile(test.boot1.2$bootstrap.adj.par[,5],probs=0.025) ## width of adjusted lembda


## Bootstrap estimate (mean) of sigma u and e for

	Bootstap.estimate.sim[i,1] =mean(test.boot1.2$bootstrap.par[,3])
	Bootstap.estimate.sim[i,2] =mean(test.boot1.2$bootstrap.par[,4])
	Bootstap.estimate.sim[i,3] =mean(test.boot1.2$bootstrap.adj.par[,3])
	Bootstap.estimate.sim[i,4] =mean(test.boot1.2$bootstrap.adj.par[,4])


	mat.Var.Sigma[i,4] <-var(test.boot1.2$bootstrap.par[,3])
	mat.Var.Sigma[i,5] <-var(test.boot1.2$bootstrap.par[,4])
	mat.Var.Sigma[i,6] <-cov(test.boot1.2$bootstrap.par[,3],test.boot1.2$bootstrap.par[,4])


	mat.Var.Sigma[i,7] <-var(test.boot1.2$bootstrap.adj.par[,3])
	mat.Var.Sigma[i,8] <-var(test.boot1.2$bootstrap.adj.par[,4])
	mat.Var.Sigma[i,9] <-cov(test.boot1.2$bootstrap.adj.par[,3],test.boot1.2$bootstrap.adj.par[,4])

## End of Simulations
#print(c("simulation number ",i))


}


rbind(round(apply(group.boot.sim1.2[sim.start:sim.finish,1:5],2,mean),4)  ## results estimates of parameters  over the simulations
,round(apply(group.boot.sim1.2[sim.start:sim.finish,1:5],2,sd),4)    ## results sd of estimates of parameters  over the simulations
,round(apply(group.boot.sim1.2[sim.start:sim.finish,6:10],2,mean)   ,2)     ## ACR for boot estimates (raw) of parameters  over the simulations
,round(apply(group.boot.sim1.2[sim.start:sim.finish,11:15],2,mean)  ,2)     ## ACR for adjusted boot estimates of parameters  over the simulations
,round(apply(group.boot.sim1.2[sim.start:sim.finish,16:20],2,mean) ,4)       ## width for boot estimates (raw) of parameters  over the simulations
,round(apply(group.boot.sim1.2[sim.start:sim.finish,21:25],2,mean),4  )     ## width for adjusted boot estimates of parameters  over the simulations
)

t(matrix(round(apply(mat.Var.Sigma[sim.start:sim.finish,],2,mean),8),3,3))
round(apply(Bootstap.estimate.sim[sim.start:sim.finish,],2,mean),4)
