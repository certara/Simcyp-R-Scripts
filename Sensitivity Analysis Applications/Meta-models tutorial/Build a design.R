library("spacefillr")

#Script to generate a design on the unit hypercube and re-scale to a model space design.
#Examples of linear and lof scaling are shown

n.sobol=2^8 #256 point design
n.param = 5 #5 parameters 

#Generate Owen scrambled Sobol sequence on [0,1]^5
Design_Unit = generate_sobol_owen_set(n.sobol, n.param)

Rescaling = data.frame(Lower = c(0.05, 0.05, 0.2, 0.55, 0.45), 
                       Upper=c(50, 10, 40, 2, 0.55), 
                       Scale=c('Log', 'Linear', 'Log', 'Linear', 'Linear'))

#Generate the model space design based on used input of lower and upper bounds and scaling.


Design_Model = matrix(rep(0,5*256), ncol=5)

#Calculate the unit design matrix based on model design points, and user design spec


for(i in 1:5)
{
  if(Rescaling$Scale[i]=="Linear")
  {
    Design_Model[,i] = Design_Unit[,i]*(Rescaling$Upper[i] - Rescaling$Lower[i]) +  Rescaling$Lower[i]
  }
  
  if(Rescaling$Scale[i]=="Log")
  {
    Design_Model[,i] = exp( Design_Unit[,i]*log(Rescaling$Upper[i] / Rescaling$Lower[i]) + log(Rescaling$Lower[i])  )
  }
  
}

#Design model is the model scale design used for model runs.
#Each row is a model scale design point corresponding to the parameters which must be changed for that run of the model.


