setwd('/Users/David/Dropbox/Oxford/Dissertation/')
getwd()

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(deSolve)
library(plot3D)

data = read.csv("datatest.csv", header=T, stringsAsFactors=F)

#Set the data to be used
inputdata = setRemissionData(data, showPlot = T)
inputdata = setPatientData(data, 25, showPlot = T)

standata = list(T = length(inputdata[,1]), 
                          R = inputdata$BCR.ABL1.ABL1....,
                          t0 = 0,
                          ts = inputdata$Month[-1])

fit <- stan(file = 'MScDissertation/odemcmc.stan', data = standata, 
            iter = 1000, chains = 1)

print(fit)
get_inits(fit)

finalsolutionplot(stanfit = fit, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
generatePlots(stanfit = fit, plotDensity = T, plotTheta = T, ploty0 = T, plotz = T, plotv = T)

svdtransfromplot(stanfit = fit, Phi = 0, Theta = 100)


#y0remissioninit = c(0.289, 0.009, 0.006, 0.04, 0.066)
#thetaremissioninit = c(0.7,0.7,0.6,0.5,0.1,0.1,0.1,0.9)

#y0relapseinit = c(0.27,0.219,0.05,0.077,0.162)
#thetarelapseinit = c(0.1,0.8,0.2,0.8,0.2,0.9,0.4,0.8)

#Full remission data, initialization for initial values
#fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_remission, iter = 1000, chains = 1,
#            init = list(list(y0temp=y0remissioninit/sum(y0remissioninit), z=sum(y0remissioninit))))
#Full remission data, full initialization
#fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_remission, iter = 100, chains = 1,
#            init = list(list(y0temp=y0remissioninit/sum(y0remissioninit), z=sum(y0remissioninit), theta=thetaremissioninit)))
#Patient 25, initial value initialization
#fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_relapse25, iter = 1000, chains = 1,
#            init = list(list(y0temp=y0relapseinit/sum(y0relapseinit), z=sum(y0relapseinit))))
#Patient 25, full initialization
#fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_relapse25, iter = 1000, chains = 1,
#            init = list(list(y0temp=y0relapseinit/sum(y0relapseinit), z=sum(y0relapseinit), theta=thetarelapseinit)))




#Write and read to file
#write.csv(extract(fit,permuted=T), file='50000relapse.csv')
mcmcdat = read.csv('50000relapse.csv', header=T)
finalsolutionplot(mcmcData = mcmcdat, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
generatePlots(mcmcData = mcmcdat, plotDensity = T, plotTheta = T, ploty0 = T, plotz = T, plotv = T)
svdtransfromplot(mcmcData = mcmcdat, Phi = 0, Theta = 100)

#Functions

#Ode system for solving ode in r
desystem = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx0 <- ax * x0 * (1 - x0 - x1 - x2 - y0 - y1) - bx * x0
    dx1 <- bx * x0 + cx * x1 * (1 - x0 - x1 - x2 - y0 - y1) - dx * x1
    dx2 <- dx * x1 - ex * x2
    dy0 <- ay * y0 * (1 - x0 - x1 - x2 - y0 - y1) - by * y0
    dy1 <- by * y0 - ey * y1
    list(c(dx0, dx1, dx2, dy0, dy1))
  })
}

#Plots the data and the ode solution using the best stan outputs
finalsolutionplot = function(stanfit = NULL, mcmcData = NULL, xdata, ydata)
{
  plot(xdata, ydata, ylim=c(0,1), xlim=c(0, max(xdata)+5))
  
  if(is.null(mcmcData))
  {
    la = extract(stanfit, permuted = T) # return a list of arrays 
    thetamcmc = la$theta
    y0mcmc = la$y0
    loglikmcmc = la$lp__
  }
  else
  {
    mcmcData = as.matrix(mcmcData)
    colnames(mcmcData) = NULL
    thetamcmc = mcmcData[,9:16]
    y0mcmc = mcmcData[,17:21]
    loglikmcmc = mcmcData[,22]
  }
  
  max_index = which.max(loglikmcmc)
  
  theta = thetamcmc[max_index,]
  y0 = y0mcmc[max_index,]
  
  print(theta)
  print(y0)
  
  parameters = c(ax = theta[1], bx = theta[2], cx = theta[3],
                 dx = theta[4], ex = theta[5], ay = theta[6],
                 by = theta[7], ey = theta[8])
  
  state = c(x0 = y0[1], x1 = y0[2], x2 = y0[3], y0 = y0[4], y1 = y0[5])
  
  times = seq(0, max(xdata) + 5, by = 0.1)
  out = ode(y = state, times = times, func = desystem, parms = parameters)
  head(out)
  
  R = out[,6]/(out[,6] + 2*out[,4])
  lines(out[,1], R, col='red')
}

#plots the svd transform of the mcmc output in 3d
svdtransfromplot = function(stanfit = NULL, mcmcData = NULL, Phi=30, Theta=30)
{
  if(is.null(mcmcData))
  {
    la = extract(stanfit, permuted = T) # return a list of arrays 
    thetamcmc = la$theta
    y0mcmc = la$y0
    vmcmc = la$v
  }
  else
  {
    mcmcData = as.matrix(mcmcData)
    colnames(mcmcData) = NULL
    thetamcmc = mcmcData[,9:16]
    y0mcmc = mcmcData[,17:21]
    vmcmc = mcmcData[,2]
  }
  
  n=14
  
  m = cbind(thetamcmc,
            y0mcmc,
            vmcmc)
  
  decomp = svd(m)
  sig = matrix(c(decomp$d[1], rep(0,n), decomp$d[2], rep(0, n), decomp$d[3], rep(0, n-3)), nrow=n)
  transform = decomp$u %*% sig
  scatter3D(transform[,1],transform[,2],transform[,3], col=1, phi=Phi, theta=Theta)
}

#Sets the patient data for a given patient
setPatientData = function(fulldata, PatientNumber = 25, showPlot = F)
{
  tempdata = fulldata
  tempdata[,3] = as.numeric(tempdata[,3])
  tempdata = na.omit(tempdata)

  tempdata[,3] = tempdata[,3]/100
  ind = which(tempdata[,1] == PatientNumber)
  
  PatientData = tempdata[ind,]
  
  if(showPlot == T)
    plot(PatientData[,2], PatientData[,3], ylim=c(0,1))
  
  return(PatientData)
}

#Sets aggrigated remission data
setRemissionData = function(fulldata, showPlot=F)
{
  tempdata = fulldata
  
  tempdata[,3] = as.numeric(tempdata[,3])
  tempdata = na.omit(tempdata)
  
  tempdata[,3] = tempdata[,3]/100
  relapse = NULL
  
  for(i in 1:(length(tempdata[,1]) - 1))
  {
    if(tempdata[i+1,3] > tempdata[i,3] + 0.05 && tempdata[i,1] == tempdata[i+1,1])
      relapse = c(relapse, tempdata[i,1])
  }
  
  relapse = unique(relapse)
  remission = setdiff(1:69,relapse)
  
  remissionDat = tempdata[which(tempdata[,1] %in% remission),]
  relapseDat = tempdata[which(tempdata[,1] %in% relapse),]
  
  remissionDat = aggregate(BCR.ABL1.ABL1.... ~ Month, data = remissionDat, mean)
  
  if(showPlot == T)
    plot(remissionDat[,1], remissionDat[,2], ylim=c(0,1))
  return(remissionDat)
}

#Generate Density and Trace plots from mcmc output
generatePlots = function(stanfit = NULL, mcmcData = NULL, plotDensity = T, plotTrace = F, plotTheta = T, ploty0 = F, plotz = F, plotv = F, plotloglik = F)
{
  if(is.null(mcmcData))
  {
    nonpermutedvalues = extract(stanfit,permuted=F)
    thetamcmc = nonpermutedvalues[,1,8:15]
    y0mcmc = nonpermutedvalues[,1,16:20]
    zmcmc = nonpermutedvalues[,1,2]
    vmcmc = nonpermutedvalues[,1,1]
    loglikmcmc = nonpermutedvalues[,1,21]
  }
  else
  {
    mcmcData = as.matrix(mcmcData)
    colnames(mcmcData) = NULL
    thetamcmc = mcmcData[,9:16]
    y0mcmc = mcmcData[,17:21]
    zmcmc = mcmcData[,3]
    vmcmc = mcmcData[,2]
    loglikmcmc = mcmcData[,22]
  }
  if(plotDensity == T)
  {
    if(plotTheta == T)
    {
      for(i in 1:8)
        plot(density(thetamcmc[,i]), main=paste("Density of theta",i))
    }
    if(ploty0 == T)
    {
      for(i in 1:5)
        plot(density(y0mcmc[,i]),main=paste("Density of y0 for ode",i))
    }
    if(plotz == T)
    {
      plot(density(zmcmc), main="Density of z")
    }
    if(plotv == T)
    {
      plot(density(vmcmc), main="Density of v")
    }
    if(plotv == T)
    {
      plot(density(loglikmcmc), main="Density of Posterior Log Likelihood")
    }
  }
  if(plotTrace == T)
  {
    if(plotTheta == T)
    {
      for(i in 1:8)
        plot(thetamcmc[,i], main=paste("Trace of theta",i), type='l')
    }
    if(ploty0 == T)
    {
      for(i in 1:5)
        plot(y0mcmc[,i],main=paste("Trace of y0 for ode",i), type='l')
    }
    if(plotz == T)
    {
      plot(zmcmc,main="Trace of z", type='l')
    }
    if(plotv == T)
    {
      plot(vmcmc,main="Trace of v", type='l')
    }
    if(plotloglik == T)
    {
      plot(loglikmcmc,main="Trace of Posterior Log Likelihood", type='l')
    }
  }
}

