

Banks = read.csv(file.path('Data_and_Model','financial_data.csv'))

source("KRSW/ModelDiagnosticsNewNames.r")

x<-Banks

# Use empirical probability integral transform to put on uniform scale

u1<-unif(x=x[,1])
hist(u1)

u2<-unif(x=x[,2])
hist(u2)

u3<-unif(x=x[,3])
hist(u3)

u4<-unif(x=x[,4])
hist(u4)

u5<-unif(x=x[,5])
hist(u5)
# Transform to standard Pareto scale

xpb<-cbind(1/(1-u1),1/(1-u2),1/(1-u3),1/(1-u4),1/(1-u5) )


# Create matrix of data on MVP scale (at least one threshold exc)
# (exponential of Y_E-u_E|Y_E \not\leq u_E)
q = 0.95
logic<-xpb[,1]>quantile(xpb[,1],q)|xpb[,2]>quantile(xpb[,2],q)|xpb[,3]>quantile(xpb[,3],q)|xpb[,4]>quantile(xpb[,4],q)|xpb[,5]>quantile(xpb[,5],q)
xpb2<-cbind(xpb[logic,1]/quantile(xpb[,1],q),xpb[logic,2]/quantile(xpb[,2],q),xpb[logic,3]/quantile(xpb[,3],q),xpb[logic,4]/quantile(xpb[,4],q),xpb[logic,5]/quantile(xpb[,5],q))
plot(log(xpb2))
dim(xpb2)

# exponential margins
xeb2<-log(xpb2)


# Same but on MGPD scale (i.e. scale of the observations)

Banks2<-cbind(Banks[logic,1]-quantile(Banks[,1],q),Banks[logic,2]-quantile(Banks[,2],q),
              Banks[logic,3]-quantile(Banks[,3],q),Banks[logic,4]-quantile(Banks[,4],q),
              Banks[logic,5]-quantile(Banks[,5],q))
plot(Banks2)
dim(Banks2)




# Initially examine "most complex" dependence models to home in on best family of models
########################################################################################
source("KRSW/CommonFunctions.r")
source("KRSW/Gumbel_T_Functions.r")
source("KRSW/MVGaussian_T_Functions.r")
source("KRSW/RevExp_T_Functions.r")
source("KRSW/RevExp_U_Functions.r")
source("KRSW/Gumbel_U_Functions.r")

#######################################################################################################
# NB: all fits below are written twice, the second version with starting values obtained from the first 
# fit. For simple models, the optimum will be obtained with a single optimization, but for others, a few
# iterations is necessary.
#######################################################################################################

# h_T based on the Gumbel

fit1<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=F, maxit=2000)
fit1
fit1<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=F, maxit=2000, dep.start=fit1$mle)
fit1

fit1.1<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=T, maxit=2000)
fit1.1
fit1.1<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=T, maxit=2000,dep.start=fit1.1$mle)
fit1.1



#####################################################################################################
#library(geoR)
library(mvtnorm)
library(Matrix)

# h_T based on the Gaussian

fit2<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,5), std=T, structured.cor=F)
# Run some repeats of this:
fit2<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,5), std=T, structured.cor=F, dep.start=fit2$mle,maxit=2000)


# Issues with "degeneracy of Nelder-Mead simplex"; BFGS gives apparent convergence
fit2<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,5), std=T, structured.cor=F, dep.start=fit2$mle, method="BFGS")

# Without free Gaussian location parameters
fit2.1<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,5), std=T, structured.cor=F, dep.start=fit2$mle,dep.loc.fix=T,maxit=2000)
fit2.1<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,5), std=T, structured.cor=F, dep.start=fit2.1$mle,
                          dep.loc.fix=T,maxit=2000,method="BFGS")


#####################################################################################################
# h_T based on the reverse exponential

fit3<-fit.MGPD.RevExpT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=F)
fit3
fit3<-fit.MGPD.RevExpT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=F, dep.start = fit3$mle)
fit3

fit3.1<-fit.MGPD.RevExpT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=T)
fit3.1
fit3.1<-fit.MGPD.RevExpT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=T, dep.start = fit3.1$mle)
fit3.1


# Minimized log-likelihoods and AICs for most complicated models

fit1$nll
fit2$nll
fit3$nll

fit1$nll+2*length(fit1$mle)
fit2$nll+2*length(fit2$mle)
fit3$nll+2*length(fit3$mle)


########################################################################################

# Gumbel model based on h_T is preferred

fit1
fit1.1

# Further parameterizations of this model

# Free scale parameter, constrained location parameters
fit1.2<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=F, dep.loc.fix=T)
fit1.2
fit1.2<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=F, dep.loc.fix=T, dep.start = fit1.2$mle)
fit1.2

# Single scale parameter, fixed location
fit1.3<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=T, dep.loc.fix=T)
fit1.3
fit1.3<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,5), std=T, dep.scale.fix=T, dep.loc.fix=T, dep.start = fit1.3$mle)
fit1.3


# NLLs to 1d.p.

# \alpha1,2,3,4; \beta1,2,3
round(fit1$nll,1)
# \alpha; \beta1,2,3
round(fit1.1$nll,1)
# \alpha1,2,3,4
round(fit1.2$nll,1)
# \alpha
round(fit1.3$nll,1)

# Likelihood ratio tests

# Sequence: 1->1.2->1.3
1-pchisq(2*(fit1.2$nll-fit1$nll), df=length(fit1$mle)-length(fit1.2$mle))
1-pchisq(2*(fit1.3$nll-fit1.2$nll), df=length(fit1.2$mle)-length(fit1.3$mle))

# Sequence 1->1.1->1.3
1-pchisq(2*(fit1.1$nll-fit1$nll), df=length(fit1$mle)-length(fit1.1$mle))
1-pchisq(2*(fit1.3$nll-fit1.1$nll), df=length(fit1.1$mle)-length(fit1.3$mle))

1-pchisq(2*(fit1.3$nll-fit1$nll), df=length(fit1$mle)-length(fit1.3$mle))



#################################################################################################
# Fit margins and dependence simulataneously
#################################################################################################

fit1.4<-fit.MGPD.GumbelT(x=Banks2, u=rep(0,5), std=F,
                         dep.scale.fix=F, dep.loc.fix=T, marg.shape.ind=1:5, marg.scale.ind=1:5,
                         dep.start=fit1.2$mle[1:5],maxit=5000,method ="BFGS" )

fit1.4<-fit.MGPD.GumbelT(x=Banks2, u=rep(0,5), std=F,
                         dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=rep(1,5), marg.scale.ind=1:5,
                         marg.scale.start=fit1.4$mle[6:10], marg.shape.start=fit1.4$mle[11:15],
                         dep.start=fit1.4$mle[1:5],maxit=5000)
fit1.4


# Test for common shape

fit1.5<-fit.MGPD.GumbelT(x=Banks2, u=rep(0,5), std=F,
                         dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=rep(1,5), marg.scale.ind=1:5,
                         marg.scale.start=fit1.4$mle[6:10], marg.shape.start=mean(fit1.4$mle[11:15]),
                         dep.start=fit1.4$mle[1:5],maxit=5000,method ="BFGS")

fit1.5

1-pchisq(2*(fit1.5$nll-fit1.4$nll),df=length(fit1.4$mle)-length(fit1.5$mle))




###############################################################################################
# Output 100 simulated datasets with estimated mGPD
###############################################################################################

install.packages("jsonlite")  # if not installed
library(jsonlite)
set.seed(1234)
par(mfrow=c(1,1))
plot(density(Banks2[,1]))

mle1 = fit1.4$mle

pred = sim.RevExpT.MGPD(10000, d= 5, a = mle1[1:5],beta = rep(0,5),
                        sig = mle1[6:10],gamma = mle1[11:15])

pred_mat = array(numeric(),c(100,100,5)) 
for (b in 1:100){
  pred = sim.RevExpT.MGPD(100, d= 5, a = mle1[1:5],beta = rep(0,5),
                          sig = mle1[6:10],gamma = mle1[11:15])
  pred_mat[b,,] = pred
}

pred_mat_json <- toJSON(pred_mat, pretty = TRUE, auto_unbox = TRUE)
writeLines(pred_mat_json, con = file.path('Data_and_Model', 'mGPD_pred.json'))

