
## This is an final script in the workflow. 
## Reads input and output files and runs meta-model based GSA

#####################################################################################################################################


Rescaling = read.csv('Theo Scales.csv')
Design = read.csv('Theo Design.csv')
Results = read.csv('Theo results.csv')


Design_Unit = matrix(rep(NA,5*256), ncol=5)

#Calculate the unit design matrix based on model design points, and user design spec
#Design_Unit


for(i in 1:5)
{
  if(Rescaling$Scale[i]=="Linear")
  {
    Design_Unit[,i] = (Design[,(i+1)] - Rescaling$Lower[i])*(Rescaling$Upper[i] - Rescaling$Lower[i])^-1
  }
  
  if(Rescaling$Scale[i]=="Log")
  {
    Design_Unit[,i] = log(Design[,(i+1)] / Rescaling$Lower[i])*log(Rescaling$Upper[i] / Rescaling$Lower[i])^-1
  }
  
}

library("sensitivity")

#This code block here generates the design points for the GSA which are on the [0,1] interval.
#We apply two transformations
#The first converts the sensitivity analysis design on the [0,1] interval to the model design space. 
#Uniform distributions are assumed throughout for model inputs. Different transformations would be applied for alternative distributions.

#The second transformation converts the model design space to the [0,1] interval for feeding into the meta model for predictions
#In the case of linear scales the transformations cancel out. For log scales the transformations do NOT cancel out.

dim=5
m = 3000

X1 <- data.frame(matrix(runif(dim * m), nrow = m))
X2 <- data.frame(matrix(runif(dim * m), nrow = m))

x = sobolSalt(model = NULL, X1, X2, scheme="A", nboot = 0)

#This contains the 6000 design points for Sobol SA
GSAdesignunit = x$X
GSAdesignmodel =  GSAdesignunit 

for(i in 1:5)

  {
  GSAdesignmodel[,i] = GSAdesignunit[,i]*(Rescaling$Upper[i] - Rescaling$Lower[i]) + Rescaling$Lower[i]
  }


SaltSobolPoints = x$X

for(i in 1:5)
{
  if(Rescaling$Scale[i]=="Linear")
  {
    SaltSobolPoints[,i] = (GSAdesignmodel[,i] - Rescaling$Lower[i])*(Rescaling$Upper[i] - Rescaling$Lower[i])^-1
  }
  
  if(Rescaling$Scale[i]=="Log")
  {
    SaltSobolPoints[,i] = log(GSAdesignmodel[,i] / Rescaling$Lower[i])*log(Rescaling$Upper[i] / Rescaling$Lower[i])^-1
  }
  
}

#We use the GSA design SaltSobolPoints at the end of GP model building 



#############################################
#############################################

#Functions for building the meta-model


corrs = function(x, b)
{
  #This is a function calculating correlations between design points
  #x is the n*p design point matrix on the unit hypercube.
  #b is the p dimensional vector of smoothing parameters.
  #Output is an n*n matrix of correlations.
  
  n = dim(x)[1]
  p = dim(x)[2]
  A = matrix(0,n,n)
  
  B = diag(b)
  
  for(i in 1:n)
  {
    for(j in 1:i)
    {
      calccorr = exp(- t( x[i,] - x[j,])  %*% B %*% (x[i,] - x[j,]))
      A[i,j] = calccorr
      A[j,i] = A[i,j]
    }
  }
  return(A)
}

regressor = function(x, func=0)
{
  #This is a function for creating the matrix of regressor variables
  #input x is an n*p matrix of points
  #func is an indicator with default 0: value of zero indicates constant mean function, value of one indicates linear mean function.
  #output is the matrix of regressor variables
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  
  if (func == 0)
  {
    out = matrix(1, n, 1)
  }
  
  if (func == 1)
  {
    out = matrix(1,n,(p+1))
    out[,2:(p+1)] = x
  } 
  
  return(out)
}


gls.estimates = function(x, Ainv, y)
{
  #This is a function for estimating statistical parameters of the GP regression model
  #input x is an n*q matrix of points where q is the length of the beta parameter vector (either length 1 or p+1)
  #y is an n*1 vector of model outputs
  #Ainv is the n*n inverse of the correlation matrix.
  #Function outputs a list with of GLS estimates for regressor coefficients beta and the variance
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  y=matrix(y,n,1)
  
  p1 = Matrix::solve(t(x) %*% Ainv %*% x)
  p2 = t(x)%*%Ainv%*%y
  
  beta = as.vector(p1%*%p2)
  
  sigmasq = 1/(n-p-2) *  t(y - x %*% beta) %*% Ainv %*% (y - x %*% beta)
  out=list(beta=beta, sigmasq=sigmasq)
  
  return(out)
  
}  


emulator.predict = function (predict, Ainv, design, scales, output, func = 0) 
  
{
  #This is a general function for making predictions at unsmampled design points
  #predict is a m*p matrix of coordinates for which predictions are required.
  #design x is an n*p matrix of points in which output is observed.
  #Ainv is the inverse of the n*n matrix of design point correlations
  #output is an n*1 vector of model outputs.
  #func is an indicator variable: zero indicates a constant mean function, value of one indicates linear mean function.
  #Function outputs a list with two elements: vector predictions; vector of prediction variances
  
  H = regressor(design, func) #returns either an n*1 or n*(p+1) matrix.
  h = regressor(predict, func) #returns either an m*1 or m*(p+1) matrix.
  
  beta = gls.estimates(H, Ainv, output)$beta
  sigmasq = gls.estimates(H, Ainv, output)$sigmasq
  
  m = dim(predict)[1]
  n = dim(design)[1]
  
  T = matrix(0,m,n)
  C = matrix(0,m,m)
  
  B = diag(scales)
  
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      calccorr =  exp(- t( predict[i,] - design[j,])  %*% B %*% (predict[i,] - design[j,]))
      T[i,j] = calccorr
    }
  }
  
  
  for(i in 1:m)
  {
    for(j in 1:i)
    {
      calccorr =  exp(- t( predict[i,] - predict[j,])  %*% B %*% (predict[i,] - predict[j,]))
      C[i,j] = calccorr
      C[j,i] = C[i,j]
    }
    
  }
  
  sigmasq = as.vector(sigmasq)
  
  predict.mean = h %*% beta + T %*% Ainv %*% (output - H %*% beta)
  
  predict.covar = sigmasq * (C - T %*% Ainv %*% t(T))
  
  out = list(means = as.vector(predict.mean),
             vars = as.vector(abs(diag(predict.covar)))
  )
  
}

cross.validate = function (input, A, output, scales, func = 0) 
{
  #This is a function for Leave one out (LOO) predictions.
  #Input is the n*p design matrix on the unit hypercube.
  #A is n*n correlation matrix
  #Output is the n*1 vector of model outputs.
  #Scales is the p dimension vector of scaling parameters
  #Func is an indicator variable: 0 corresponds to constant mean function, 1 corresponds to a linear mean function.
  #Function outputs a list with two elements: vector of LOO predictions; vector of prediction variances
  
  d.missing.estimated <- output + NA
  d.missing.var <- output + NA
  for (i in 1:nrow(input)) {
    val.oneshort <- input[-i, , drop = FALSE]
    val.missing <- as.matrix(t(input[i, ]))
    d.oneshort <- output[-i]
    d.missing <- output[i]
    A.oneshort <- A[-i, -i]
    Ainv.oneshort <- Matrix::solve(A.oneshort)
    
    d.missing.estimated[i]  <- emulator.predict(val.missing, Ainv.oneshort, val.oneshort, scales, d.oneshort, func)$means
    d.missing.var[i] <- emulator.predict(val.missing, Ainv.oneshort, val.oneshort, scales, d.oneshort, func)$vars
  }
  
  out = list(means=d.missing.estimated, 
             vars=d.missing.var)
  
  return(out)
}

objectivefunc_constant <- function(p)
{
  #This function for optimization assumes a mean function that is constant
  #Linear mean is specified through the func = 0 specification in the cross.validate function
  #Sole input is the vector of scale parameters for the correlation function.
  #Function output is the (negative) sum of squared LOO residuals
  smooth = exp(p)
  A = corrs(x=model.input, smooth)
  n = nrow(A)
  A = A + 0.00000001*diag(n)
  crossval = cross.validate(model.input, A, model.output, scales=smooth, func=0)$means
  error = crossval - model.output
  out = sum(error^2)
  return(out)
}


objectivefunc_constant <- function(p)
{
  #This function for optimization assumes a mean function that is constant
  #Linear mean is specified through the func = 0 specification in the cross.validate function
  #Sole input is the vector of scale parameters for the correlation function.
  #Function output is the (negative) sum of squared LOO residuals
  smooth = exp(p)
  A = corrs(x=model.input, smooth)
  n = nrow(A)
  q=1
  A = A + 0.00000001*diag(n)
  Ainv = solve(A)
  H = matrix(1, n, 1)
  beta = gls.estimates(H, Ainv, model.output)$beta
  sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq
  
  p1 = - 0.5*(n-q)*log(sqrt(sigmasq)) 
  p2 = - 0.5*log(det(A))
  p3 = - 0.5*log(det(t(H)%*%A%*%H))
  
  out = (p1 + p2 + p3)/-1
  return(out)
}



emulator.predict.fast = function (predict, Ainv, design, scales, output, func = 0) 
  
{
  #This is a general function for making predictions at unsmampled design points
  #predict is a m*p matrix of coordinates for which predictions are required.
  #design x is an n*p matrix of points in which output is observed.
  #Ainv is the inverse of the n*n matrix of design point correlations
  #output is an n*1 vector of model outputs.
  #func is an indicator variable: zero indicates a constant mean function, value of one indicates linear mean function.
  #Function outputs a list with two elements: vector predictions; vector of prediction variances
  
  H = regressor(design, func) #returns either an n*1 or n*(p+1) matrix.
  h = regressor(predict, func) #returns either an m*1 or m*(p+1) matrix.
  
  beta = gls.estimates(H, Ainv, output)$beta
  sigmasq = gls.estimates(H, Ainv, output)$sigmasq
  
  m = dim(predict)[1]
  n = dim(design)[1]
  
  T = matrix(0,m,n)
  C = matrix(0,m,m)
  
  B = diag(scales)
  
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      calccorr =  exp(- t( predict[i,] - design[j,])  %*% B %*% (predict[i,] - design[j,]))
      T[i,j] = calccorr
    }
  }
  
  predict.mean = h %*% beta + T %*% Ainv %*% (output - H %*% beta)
  
  out = list(means = as.vector(predict.mean))
  
  
}

#####################################################
#####################################################


model.output <- as.vector(Results$AUC[1:64])
model.input <- Design_Unit[1:64,]


nparam=5

initial.values <- runif(5, 0, 1)

g <- optim(log(initial.values),objectivefunc_constant,method="SANN" ,control=list(trace=0,maxit=500,fnscale= 1))
g


optim.scales = exp(g$par)

A <- corrs(x=model.input, optim.scales)
n = nrow(A)
A = A + 0.00000001*diag(n)
crossval.m = cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$means
crossval.sd = sqrt(cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$vars)
error = crossval.m - model.output
out = - sum(error^2)
out

plot_AUC_results = data.frame(model.output, crossval.m, 
                              crossval.lower = crossval.m - 2*crossval.sd, crossval.upper = crossval.m + 2*crossval.sd, error)

write.csv(plot_AUC_results, 'AUC linear_64.csv')

optim.scales.AUC = optim.scales
A.AUC = A

Ainv = Matrix::solve(A.AUC)
H = regressor(model.input, func=0)
beta = gls.estimates(H, Ainv, model.output)$beta
sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq


mu = mean(model.output)
SS_AUC = sum((model.output - crossval.m)^2)
Q2_AUC = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_AUC

library(ggplot2)

plotres = ggplot(data=plot_AUC_results, aes(x=model.output, y=crossval.m)) +
  geom_errorbar(aes(ymin=crossval.lower, ymax=crossval.upper), width=.1) +
  geom_line(aes(x=model.output, y=model.output)) +
  geom_point(aes(x=model.output, y=crossval.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'AUC, n=64, Q2 = 0.996')

ggsave(plotres, filename= 'AUC linear_64.png')


mu = mean(model.output)
SS_AUC = sum((model.output - crossval.m)^2)
Q2_AUC = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_AUC


model.pred <- as.vector(Results$AUC[65:128])
samples1 <- Design_Unit[65:128,]

scales = optim.scales.AUC
Ainv = Matrix::solve(A.AUC)

out.of.sample.m = emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$means
out.of.sample.sd = sqrt(emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$vars)


plot_AUC_oos_constant = data.frame(model.pred, out.of.sample.m, 
                                   out.of.sample.m.lower = out.of.sample.m - 2*out.of.sample.sd, 
                                   out.of.sample.m.upper = out.of.sample.m + 2*out.of.sample.sd)

write.csv(plot_AUC_oos_constant, 'AUC OOS predictions linear_64.csv')

mu = mean(model.pred)
SS_AUC = sum((model.pred - out.of.sample.m)^2)
Q2_AUC = 1 - sum((model.pred - out.of.sample.m)^2)/sum((model.pred - mu)^2)
Q2_AUC

plotresAUC64 = ggplot(data=plot_AUC_oos_constant, aes(x=model.pred, y=out.of.sample.m)) +
  geom_errorbar(aes(ymin=out.of.sample.m.lower, ymax=out.of.sample.m.upper), width=.1) +
  geom_line(aes(x=model.pred, y=model.pred)) +
  geom_point(aes(x=model.pred, y=out.of.sample.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  scale_x_continuous(limits=c(0,70)) +
  scale_y_continuous(limits=c(0,70)) +
  labs(subtitle = 'AUC, n=64, Q2 = 0.985')

ggsave(plotres, filename='AUC OOS predictions linear_64.png')



###########################################
###########################################

#Refine

model.output <- as.vector(Results$AUC[1:128])
model.input <- Design_Unit[1:128,]


nparam=5

initial.values <- runif(5, 0, 1)

g <- optim(log(initial.values),objectivefunc_constant,method="SANN" ,control=list(trace=0,maxit=500,fnscale= 1))
g

optim.scales = exp(g$par)


A <- corrs(x=model.input, optim.scales)
n = nrow(A)
A = A + 0.00000001*diag(n)
crossval.m = cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$means
crossval.sd = sqrt(cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$vars)
error = crossval.m - model.output
out = - sum(error^2)
out



plot_AUC_results = data.frame(model.output, crossval.m, 
                              crossval.lower = crossval.m - 2*crossval.sd, crossval.upper = crossval.m + 2*crossval.sd, error)

write.csv(plot_AUC_results, 'AUC linear_128.csv')

optim.scales.AUC = optim.scales
A.AUC = A

Ainv = Matrix::solve(A)
H = regressor(model.input, func=0)
beta = gls.estimates(H, Ainv, model.output)$beta
sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq


mu = mean(model.output)
SS_AUC = sum((model.output - crossval.m)^2)
Q2_AUC = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_AUC

plotres = ggplot(data=plot_AUC_results, aes(x=model.output, y=crossval.m)) +
  geom_errorbar(aes(ymin=crossval.lower, ymax=crossval.upper), width=.1) +
  geom_line(aes(x=model.output, y=model.output)) +
  geom_point(aes(x=model.output, y=crossval.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'AUC, n=128, Q2 = 0.992')

ggsave(plotres, filename= 'AUC linear_128.png')




model.pred <- as.vector(Results$AUC[129:256])
samples1 <- Design_Unit[129:256,]

scales = optim.scales.AUC
Ainv = Matrix::solve(A.AUC)

out.of.sample.m = emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$means
out.of.sample.sd = sqrt(emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$vars)


plot_AUC_oos_constant = data.frame(model.pred, out.of.sample.m, 
                                   out.of.sample.m.lower = out.of.sample.m - 2*out.of.sample.sd, 
                                   out.of.sample.m.upper = out.of.sample.m + 2*out.of.sample.sd)


write.csv(plot_AUC_oos_constant, 'AUC OOS predictions linear_128.csv')

mu = mean(model.pred)
SS_AUC = sum((model.pred - out.of.sample.m)^2)
Q2_AUC = 1 - sum((model.pred - out.of.sample.m)^2)/sum((model.pred - mu)^2)
Q2_AUC

plotres = ggplot(data=plot_AUC_oos_constant, aes(x=model.pred, y=out.of.sample.m)) +
  geom_errorbar(aes(ymin=out.of.sample.m.lower, ymax=out.of.sample.m.upper), width=.1) +
  geom_line(aes(x=model.pred, y=model.pred)) +
  geom_point(aes(x=model.pred, y=out.of.sample.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'AUC, n=128, Q2 = 0.994')

ggsave(plotres, filename='AUC OOS predictions linear_128.png')




######################################################################
######################################################################

#Predictions at the GSA design points

predict = emulator.predict.fast(SaltSobolPoints, Ainv, model.input, optim.scales, model.output, func = 0)$means

#Pass to the Sobol GSA object and gather

tell(x, predict)
x$V
main = x$S
total = x$T

AUCres128 = data.frame(main, total)
write.csv(AUCres128, 'AUCres128.csv')


#########################################################################
#########################################################################

#Refine

model.output <- as.vector(Results$AUC[1:256])
model.input <- Design_Unit[1:256,]


nparam=5
initial.values <- runif(5, 0, 1)

g <- optim(log(initial.values),objectivefunc_constant,method="SANN" ,control=list(trace=0,maxit=500,fnscale= 1))
g

optim.scales = exp(g$par)


A <- corrs(x=model.input, optim.scales)
n = nrow(A)
A = A + 0.00000001*diag(n)
crossval.m = cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$means
crossval.sd = sqrt(cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$vars)
error = crossval.m - model.output
out = - sum(error^2)
out



plot_AUC_results = data.frame(model.output, crossval.m, 
                              crossval.lower = crossval.m - 2*crossval.sd, crossval.upper = crossval.m + 2*crossval.sd, error)

write.csv(plot_AUC_results, 'AUC linear_256.csv')

optim.scales.AUC = optim.scales
A.AUC = A

Ainv = Matrix::solve(A)
H = regressor(model.input, func=0)
beta = gls.estimates(H, Ainv, model.output)$beta
sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq

plotres = ggplot(data=plot_AUC_results, aes(x=model.output, y=crossval.m)) +
  geom_errorbar(aes(ymin=crossval.lower, ymax=crossval.upper), width=.1) +
  geom_line(aes(x=model.output, y=model.output)) +
  geom_point(aes(x=model.output, y=crossval.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'AUC, n=128, Q2 = 0.993')

ggsave(plotres, filename= 'AUC linear_256.png')


mu = mean(model.output)
SS_AUC = sum((model.output - crossval.m)^2)
Q2_AUC = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_AUC

######################################################################
######################################################################

#Predictions at the GSA design points

predict = emulator.predict.fast(SaltSobolPoints, Ainv, model.input, optim.scales, model.output, func = 0)$means

#Pass to the Sobol object and gather

tell(x, predict)
x$V
main = x$S
total = x$T

AUCres256 = data.frame(main, total)
write.csv(AUCres256, 'AUCres256.csv')


#########################################################################
#########################################################################


model.output <- as.vector(Results$Cmax[1:64])
model.input <- Design_Unit[1:64,]


nparam=5
initial.values <- runif(5, 0, 1)

g <- optim(log(initial.values),objectivefunc_constant,method="SANN" ,control=list(trace=0,maxit=1500,fnscale= 1))
g


optim.scales = exp(g$par)

A <- corrs(x=model.input, optim.scales)
n = nrow(A)
A = A + 0.00000001*diag(n)
crossval.m = cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$means
crossval.sd = sqrt(cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$vars)
error = crossval.m - model.output
out = - sum(error^2)
out

plot_Cmax_results = data.frame(model.output, crossval.m, 
                               crossval.lower = crossval.m - 2*crossval.sd, crossval.upper = crossval.m + 2*crossval.sd, error)

write.csv(plot_Cmax_results, 'Cmax linear_64.csv')

optim.scales.Cmax = optim.scales
A.Cmax = A

Ainv = Matrix::solve(A.Cmax)
H = regressor(model.input, func=0)
beta = gls.estimates(H, Ainv, model.output)$beta
sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq


mu = mean(model.output)
SS_Cmax = sum((model.output - crossval.m)^2)
Q2_Cmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Cmax

library(ggplot2)

plotres = ggplot(data=plot_Cmax_results, aes(x=model.output, y=crossval.m)) +
  geom_errorbar(aes(ymin=crossval.lower, ymax=crossval.upper), width=.1) +
  geom_line(aes(x=model.output, y=model.output)) +
  geom_point(aes(x=model.output, y=crossval.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Cmax, n=64, Q2 = 0.987')

ggsave(plotres, filename= 'Cmax linear_64.png')




model.pred <- as.vector(Results$Cmax[65:128])
samples1 <- Design_Unit[65:128,]

scales = optim.scales.Cmax
Ainv = Matrix::solve(A.Cmax)

out.of.sample.m = emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$means
out.of.sample.sd = sqrt(emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$vars)


plot_Cmax_oos_constant = data.frame(model.pred, out.of.sample.m, 
                                    out.of.sample.m.lower = out.of.sample.m - 2*out.of.sample.sd, 
                                    out.of.sample.m.upper = out.of.sample.m + 2*out.of.sample.sd)

write.csv(plot_Cmax_oos_constant, 'Cmax OOS predictions linear_64.csv')

mu = mean(model.pred)
SS_Cmax = sum((model.pred - out.of.sample.m)^2)
Q2_Cmax = 1 - sum((model.pred - out.of.sample.m)^2)/sum((model.pred - mu)^2)
Q2_Cmax

plotres = ggplot(data=plot_Cmax_oos_constant, aes(x=model.pred, y=out.of.sample.m)) +
  geom_errorbar(aes(ymin=out.of.sample.m.lower, ymax=out.of.sample.m.upper), width=.1) +
  geom_line(aes(x=model.pred, y=model.pred)) +
  geom_point(aes(x=model.pred, y=out.of.sample.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Cmax, n=64, Q2 = 0.949')

ggsave(plotres, filename='Cmax OOS predictions linear_64.png')


###########################################
###########################################

#Refine

model.output <- as.vector(Results$Cmax[1:128])
model.input <- Design_Unit[1:128,]


nparam=5
initial.values <- runif(5, 0, 1)

g <- optim(log(initial.values),objectivefunc_constant,method="SANN" ,control=list(trace=0,maxit=1500,fnscale= 1))
g

optim.scales = exp(g$par)


A <- corrs(x=model.input, optim.scales)
n = nrow(A)
A = A + 0.00000001*diag(n)
crossval.m = cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$means
crossval.sd = sqrt(cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$vars)
error = crossval.m - model.output
out = - sum(error^2)
out

mu = mean(model.output)
SS_Cmax = sum((model.output - crossval.m)^2)
Q2_Cmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Cmax

plot_Cmax_results = data.frame(model.output, crossval.m, 
                               crossval.lower = crossval.m - 2*crossval.sd, crossval.upper = crossval.m + 2*crossval.sd, error)

write.csv(plot_Cmax_results, 'Cmax linear_128.csv')

optim.scales.Cmax = optim.scales
A.Cmax = A

Ainv = Matrix::solve(A)
H = regressor(model.input, func=0)
beta = gls.estimates(H, Ainv, model.output)$beta
sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq

mu = mean(model.output)
SS_Cmax = sum((model.output - crossval.m)^2)
Q2_Cmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Cmax

plotres = ggplot(data=plot_Cmax_results, aes(x=model.output, y=crossval.m)) +
  geom_errorbar(aes(ymin=crossval.lower, ymax=crossval.upper), width=.1) +
  geom_line(aes(x=model.output, y=model.output)) +
  geom_point(aes(x=model.output, y=crossval.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Cmax, n=128, Q2 = 0.973')

ggsave(plotres, filename= 'Cmax linear_128.png')


model.pred <- as.vector(Results$Cmax[129:256])
samples1 <- Design_Unit[129:256,]

scales = optim.scales.Cmax
Ainv = Matrix::solve(A.Cmax)

out.of.sample.m = emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$means
out.of.sample.sd = sqrt(emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$vars)


plot_Cmax_oos_constant = data.frame(model.pred, out.of.sample.m, 
                                    out.of.sample.m.lower = out.of.sample.m - 2*out.of.sample.sd, 
                                    out.of.sample.m.upper = out.of.sample.m + 2*out.of.sample.sd)


write.csv(plot_Cmax_oos_constant, 'Cmax OOS predictions linear_128.csv')

mu = mean(model.pred)
SS_Cmax = sum((model.pred - out.of.sample.m)^2)
Q2_Cmax = 1 - sum((model.pred - out.of.sample.m)^2)/sum((model.pred - mu)^2)
Q2_Cmax

plotres = ggplot(data=plot_Cmax_oos_constant, aes(x=model.pred, y=out.of.sample.m)) +
  geom_errorbar(aes(ymin=out.of.sample.m.lower, ymax=out.of.sample.m.upper), width=.1) +
  geom_line(aes(x=model.pred, y=model.pred)) +
  geom_point(aes(x=model.pred, y=out.of.sample.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Cmax, n=128, Q2 = 0.935')

ggsave(plotres, filename='Cmax OOS predictions linear_128.png')



######################################################################
######################################################################

#Predictions at the GSA design points

predict = emulator.predict.fast(SaltSobolPoints, Ainv, model.input, optim.scales, model.output, func = 0)$means

#Pass to the Sobol GSA object and gather

tell(x, predict)
x$V
main = x$S
total = x$T

Cmaxres128 = data.frame(main, total)
write.csv(Cmaxres128, 'Cmaxres128.csv')


#########################################################################
#########################################################################

#Refine

model.output <- as.vector(Results$Cmax[1:256])
model.input <- Design_Unit[1:256,]


nparam=5
#initial.values <- c(2.510746e+00, 2.085357e+01, 1.717398e-02, 9.723053e-04, 1.591372e-02)
initial.values <- runif(5, 0, 5)

g <- optim(log(initial.values),objectivefunc_constant,method="Nelder-Mead" ,control=list(trace=0,maxit=500,fnscale= 1))
g

optim.scales = exp(g$par)


A <- corrs(x=model.input, optim.scales)
n = nrow(A)
A = A + 0.00000001*diag(n)
crossval.m = cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$means
crossval.sd = sqrt(cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$vars)
error = crossval.m - model.output
out = - sum(error^2)
out


mu = mean(model.output)
SS_Cmax = sum((model.output - crossval.m)^2)
Q2_Cmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Cmax


plot_Cmax_results = data.frame(model.output, crossval.m, 
                               crossval.lower = crossval.m - 2*crossval.sd, crossval.upper = crossval.m + 2*crossval.sd, error)

write.csv(plot_Cmax_results, 'Cmax linear_256.csv')

optim.scales.Cmax = optim.scales
A.Cmax = A

Ainv = Matrix::solve(A)
H = regressor(model.input, func=0)
beta = gls.estimates(H, Ainv, model.output)$beta
sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq

plotres = ggplot(data=plot_Cmax_results, aes(x=model.output, y=crossval.m)) +
  geom_errorbar(aes(ymin=crossval.lower, ymax=crossval.upper), width=.1) +
  geom_line(aes(x=model.output, y=model.output)) +
  geom_point(aes(x=model.output, y=crossval.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Cmax, n=128, Q2 = 0.957')

ggsave(plotres, filename= 'Cmax linear_256.png')


mu = mean(model.output)
SS_Cmax = sum((model.output - crossval.m)^2)
Q2_Cmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Cmax

######################################################################
######################################################################

#Predictions at the GSA design points

predict = emulator.predict.fast(SaltSobolPoints, Ainv, model.input, optim.scales, model.output, func = 0)$means

#Pass to the Sobol object and gather

tell(x, predict)
x$V
main = x$S
total = x$T

Cmaxres256 = data.frame(main, total)
write.csv(Cmaxres256, 'Cmaxres256.csv')


#########################################################################
#########################################################################


model.output <- as.vector(Results$Tmax[1:64])
model.input <- Design_Unit[1:64,]


nparam=5
initial.values <- runif(5, 0, 1)

g <- optim(log(initial.values),objectivefunc_constant,method="SANN" ,control=list(trace=0,maxit=1500,fnscale= 1))
g


optim.scales = exp(g$par)

A <- corrs(x=model.input, optim.scales)
n = nrow(A)
A = A + 0.00000001*diag(n)
crossval.m = cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$means
crossval.sd = sqrt(cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$vars)
error = crossval.m - model.output
out = - sum(error^2)
out

plot_Tmax_results = data.frame(model.output, crossval.m, 
                               crossval.lower = crossval.m - 2*crossval.sd, crossval.upper = crossval.m + 2*crossval.sd, error)

write.csv(plot_Tmax_results, 'Tmax linear_64.csv')

optim.scales.Tmax = optim.scales
A.Tmax = A

Ainv = Matrix::solve(A.Tmax)
H = regressor(model.input, func=0)
beta = gls.estimates(H, Ainv, model.output)$beta
sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq


mu = mean(model.output)
SS_Tmax = sum((model.output - crossval.m)^2)
Q2_Tmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Tmax

library(ggplot2)

plotres = ggplot(data=plot_Tmax_results, aes(x=model.output, y=crossval.m)) +
  geom_errorbar(aes(ymin=crossval.lower, ymax=crossval.upper), width=.1) +
  geom_line(aes(x=model.output, y=model.output)) +
  geom_point(aes(x=model.output, y=crossval.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Tmax, n=64, Q2 = 0.987')

ggsave(plotres, filename= 'Tmax linear_64.png')




model.pred <- as.vector(Results$Tmax[65:128])
samples1 <- Design_Unit[65:128,]

scales = optim.scales.Tmax
Ainv = Matrix::solve(A.Tmax)

out.of.sample.m = emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$means
out.of.sample.sd = sqrt(emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$vars)


plot_Tmax_oos_constant = data.frame(model.pred, out.of.sample.m, 
                                    out.of.sample.m.lower = out.of.sample.m - 2*out.of.sample.sd, 
                                    out.of.sample.m.upper = out.of.sample.m + 2*out.of.sample.sd)

write.csv(plot_Tmax_oos_constant, 'Tmax OOS predictions linear_64.csv')

mu = mean(model.pred)
SS_Tmax = sum((model.pred - out.of.sample.m)^2)
Q2_Tmax = 1 - sum((model.pred - out.of.sample.m)^2)/sum((model.pred - mu)^2)
Q2_Tmax

plotres = ggplot(data=plot_Tmax_oos_constant, aes(x=model.pred, y=out.of.sample.m)) +
  geom_errorbar(aes(ymin=out.of.sample.m.lower, ymax=out.of.sample.m.upper), width=.1) +
  geom_line(aes(x=model.pred, y=model.pred)) +
  geom_point(aes(x=model.pred, y=out.of.sample.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Tmax, n=64, Q2 = 0.949')

ggsave(plotres, filename='Tmax OOS predictions linear_64.png')


###########################################
###########################################

#Refine

model.output <- as.vector(Results$Tmax[1:128])
model.input <- Design_Unit[1:128,]


nparam=5
initial.values <- runif(5, 0, 1)

g <- optim(log(initial.values),objectivefunc_constant,method="SANN" ,control=list(trace=0,maxit=1500,fnscale= 1))
g

optim.scales = exp(g$par)


A <- corrs(x=model.input, optim.scales)
n = nrow(A)
A = A + 0.00000001*diag(n)
crossval.m = cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$means
crossval.sd = sqrt(cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$vars)
error = crossval.m - model.output
out = - sum(error^2)
out

mu = mean(model.output)
SS_Tmax = sum((model.output - crossval.m)^2)
Q2_Tmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Tmax

plot_Tmax_results = data.frame(model.output, crossval.m, 
                               crossval.lower = crossval.m - 2*crossval.sd, crossval.upper = crossval.m + 2*crossval.sd, error)

write.csv(plot_Tmax_results, 'Tmax linear_128.csv')

optim.scales.Tmax = optim.scales
A.Tmax = A

Ainv = Matrix::solve(A)
H = regressor(model.input, func=0)
beta = gls.estimates(H, Ainv, model.output)$beta
sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq

mu = mean(model.output)
SS_Tmax = sum((model.output - crossval.m)^2)
Q2_Tmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Tmax

plotres = ggplot(data=plot_Tmax_results, aes(x=model.output, y=crossval.m)) +
  geom_errorbar(aes(ymin=crossval.lower, ymax=crossval.upper), width=.1) +
  geom_line(aes(x=model.output, y=model.output)) +
  geom_point(aes(x=model.output, y=crossval.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Tmax, n=128, Q2 = 0.973')

ggsave(plotres, filename= 'Tmax linear_128.png')


model.pred <- as.vector(Results$Tmax[129:256])
samples1 <- Design_Unit[129:256,]

scales = optim.scales.Tmax
Ainv = Matrix::solve(A.Tmax)

out.of.sample.m = emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$means
out.of.sample.sd = sqrt(emulator.predict(samples1, Ainv, model.input, scales, model.output, func = 0)$vars)


plot_Tmax_oos_constant = data.frame(model.pred, out.of.sample.m, 
                                    out.of.sample.m.lower = out.of.sample.m - 2*out.of.sample.sd, 
                                    out.of.sample.m.upper = out.of.sample.m + 2*out.of.sample.sd)


write.csv(plot_Tmax_oos_constant, 'Tmax OOS predictions linear_128.csv')

mu = mean(model.pred)
SS_Tmax = sum((model.pred - out.of.sample.m)^2)
Q2_Tmax = 1 - sum((model.pred - out.of.sample.m)^2)/sum((model.pred - mu)^2)
Q2_Tmax

plotres = ggplot(data=plot_Tmax_oos_constant, aes(x=model.pred, y=out.of.sample.m)) +
  geom_errorbar(aes(ymin=out.of.sample.m.lower, ymax=out.of.sample.m.upper), width=.1) +
  geom_line(aes(x=model.pred, y=model.pred)) +
  geom_point(aes(x=model.pred, y=out.of.sample.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Tmax, n=128, Q2 = 0.935')

ggsave(plotres, filename='Tmax OOS predictions linear_128.png')



######################################################################
######################################################################

#Predictions at the GSA design points

predict = emulator.predict.fast(SaltSobolPoints, Ainv, model.input, optim.scales, model.output, func = 0)$means

#Pass to the Sobol GSA object and gather

tell(x, predict)
x$V
main = x$S
total = x$T

Tmaxres128 = data.frame(main, total)
write.csv(Tmaxres128, 'Tmaxres128.csv')


#########################################################################
#########################################################################

#Refine

model.output <- as.vector(Results$Tmax[1:256])
model.input <- Design_Unit[1:256,]


nparam=5
#initial.values <- c(2.510746e+00, 2.085357e+01, 1.717398e-02, 9.723053e-04, 1.591372e-02)
initial.values <- runif(5, 0, 5)

g <- optim(log(initial.values),objectivefunc_constant,method="Nelder-Mead" ,control=list(trace=0,maxit=500,fnscale= 1))
g

optim.scales = exp(g$par)


A <- corrs(x=model.input, optim.scales)
n = nrow(A)
A = A + 0.00000001*diag(n)
crossval.m = cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$means
crossval.sd = sqrt(cross.validate(model.input, A, model.output, scales=optim.scales, func=0)$vars)
error = crossval.m - model.output
out = - sum(error^2)
out


mu = mean(model.output)
SS_Tmax = sum((model.output - crossval.m)^2)
Q2_Tmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Tmax


plot_Tmax_results = data.frame(model.output, crossval.m, 
                               crossval.lower = crossval.m - 2*crossval.sd, crossval.upper = crossval.m + 2*crossval.sd, error)

write.csv(plot_Tmax_results, 'Tmax linear_256.csv')

optim.scales.Tmax = optim.scales
A.Tmax = A

Ainv = Matrix::solve(A)
H = regressor(model.input, func=0)
beta = gls.estimates(H, Ainv, model.output)$beta
sigmasq = gls.estimates(H, Ainv, model.output)$sigmasq

plotres = ggplot(data=plot_Tmax_results, aes(x=model.output, y=crossval.m)) +
  geom_errorbar(aes(ymin=crossval.lower, ymax=crossval.upper), width=.1) +
  geom_line(aes(x=model.output, y=model.output)) +
  geom_point(aes(x=model.output, y=crossval.m)) +
  xlab('Model output') +
  ylab('Emulator prediction') +
  labs(subtitle = 'Tmax, n=128, Q2 = 0.957')

ggsave(plotres, filename= 'Tmax linear_256.png')


mu = mean(model.output)
SS_Tmax = sum((model.output - crossval.m)^2)
Q2_Tmax = 1 - sum((model.output - crossval.m)^2)/sum((model.output - mu)^2)
Q2_Tmax

######################################################################
######################################################################

#Predictions at the GSA design points

predict = emulator.predict.fast(SaltSobolPoints, Ainv, model.input, optim.scales, model.output, func = 0)$means

#Pass to the Sobol object and gather

tell(x, predict)
x$V
main = x$S
total = x$T

Tmaxres256 = data.frame(main, total)
write.csv(Tmaxres256, 'Tmaxres256.csv')


