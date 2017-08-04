setwd('/Users/David/Dropbox/Oxford/Dissertation/')

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(deSolve)
library(plot3D)
library(rBeta2009)
library(scatterplot3d)
library(shinystan)

data = read.csv("datatest.csv", header=T, stringsAsFactors=F)

#Set the data to be used
inputdata = setRemissionData(data, showPlot = T)

inputdata = setPatientData(data, 10, showPlot = F)
inputdata = setPatientData(data, 25, showPlot = F)

standata = list(T = length(inputdata[,1]), 
                R = inputdata$BCR.ABL1.ABL1....,
                t0 = 0,
                ts = inputdata$Month[-1],
                alpha = c(0.9,0.9,0.4,0.9,1.5))

#number of mcmc chains
n = 4

fit <- stan(file = 'MScDissertation/odemcmc.stan', data = standata, 
              iter = 30000, chains = n, warmup = 500)

#sso = as.shinystan(samplefit)
#launch_shinystan(sso)

#Create list of mcmc output from 'fit' for use in functions
sampleMcmcList = generateListData(fit, nruns = n)
sampleMcmcList = generateListData(samplefit, nruns = n)


par(mfrow = c(1,2))
finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1...., meanSol = T)
generatePlots(mcmcList = sampleMcmcList, plotDensity = T, plotTheta = T, ploty0 = T, plotz = F, plotv = T, plotR0 = T, plotloglik = T)
generatePlots(mcmcList = sampleMcmcList, plotDensity = F, plotTrace = T, plotTheta = F, ploty0 = F, plotz = F, plotv = F, plotloglik = T)
IndividualSolutionPlot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
PlotSampleSolutions(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1...., numSamples = 30, plotIndSamples = F)
svdtransfromplot(mcmcList  = sampleMcmcList, angle = 20, dim=2)
generatePlots(mcmcList = sampleMcmcList, plotDensity = F, plotTrace = T, plotTheta = T, ploty0 = T, plotz = F, plotv = T, plotloglik = T)
generatediagnostics(mcmcfit = samplefit, plotTrace = T, plotACF = T, plotDensity = F, plotTheta = F, ploty0 = F, plotv = F, plotloglik = T)


#Running multiple patients
PatientIDs = c(4,5)

for(i in PatientIDs)
{
  inputdata = setPatientData(data, i, showPlot = F)
  
  standata = list(T = length(inputdata[,1]), 
                  R = inputdata$BCR.ABL1.ABL1....,
                  t0 = 0,
                  ts = inputdata$Month[-1],
                  alpha = c(0.9,0.9,0.4,0.9,1.5))

  n = 4
  print(paste("##########################        PATIENT",i, "     ##########################"), quote = F)
  
  fit <- stan(file = 'MScDissertation/odemcmc.stan', data = standata, 
              iter = 10000, chains = n, warmup = 500)
  
  saveRDS(fit, file = paste(i, "10000fit.RData", sep="_"))
}

for(i in PatientIDs)
{
  samplefit = readRDS(paste(i, "10000fit.RData", sep="_"))
  
  inputdata = setPatientData(data, i, showPlot = F)
  
  sampleMcmcList = generateListData(samplefit, nruns = n)
  
  finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
  generatediagnostics(mcmcfit = samplefit, plotTrace = F, plotACF = T, plotDensity = F, plotTheta = F, ploty0 = F, plotv = F, plotloglik = T)
}

#Prior Elic
priorSampling(n = 50, plotTraces = T, plotDensities = F, alpha = c(0.9,0.9,0.4,0.9,1.5))
priorSampling(n = 1000000, plotTraces = F, plotDensities = T, alpha = c(0.9,0.9,0.4,0.9,1.5)) 
  
sigsamp = rgamma(1000000,1,2)
rsamp = runif(1000000,0,1)
v = (exp(sigsamp^2) - 1)*rsamp^2
sum(v>100)
v = v[-which(v > 10)]
median(sqrt(v))

plot(density(sqrt(v)), xlim = c(0,1), main="density of sampled standard deviation", xlab = "sqrt(v)")
  

#EDA

remmissionInd = c(10, 23, 48)
relapseInd = c(7, 25, 47)

par(mfrow = c(2,3))
dataPlot(data, PatientInd = c(remmissionInd, relapseInd), overlay = F, type = 'l')

#Write and read to file
saveRDS(fit, file = "10Remission30000Normalfit.RData")
samplefit = readRDS("10Remission30000Normalfit.RData")
sampleMcmcList = readRDS("25Relapse30000_3000Normal.RData")
sampleMcmcList = readRDS("10Remission30000_3000Normal.RData")






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


generateListData = function(stanfit = NULL, mcmcData = NULL, nruns = 1)
{
  tempList = vector("list", nruns)
  
  if(is.null(mcmcData))
  {
    la = as.array(stanfit)

    for(i in seq_len(nruns))
      tempList[[i]] = la[,i,]
  }
  else
  {
    tempList[[1]] = mcmcData
  }
  return(tempList)
}

generatediagnostics = function(mcmcfit = NULL, plotTrace = F, plotACF = F, plotDensity = F, plotTheta = F, ploty0 = F, plotv = F, plotloglik = F)
{
  print(mcmcfit, pars = c("theta", "y0", "v", "lp__"), probs=NULL)
  
  thetaLabel = NULL
  y0Label = NULL
  vLabel = NULL
  logpostLabel = NULL
  
  if(plotTheta == T)
    thetaLabel = "theta"
  if(ploty0 == T)
    y0Label = "y0"
  if(plotv == T)
    vLabel = "v"
  if(plotloglik == T)
    logpostLabel = "lp__"
  
  if(plotTrace == T)
    print(stan_trace(mcmcfit, pars = c(thetaLabel, y0Label, vLabel, logpostLabel), inc_warmup = TRUE, nrow = 2) + xlim(0, 5000)) 
  
  if(plotACF == T)
    print(stan_ac(mcmcfit, pars = c(thetaLabel, y0Label, vLabel, logpostLabel)))
  
  if(plotDensity == T)
    print(stan_dens(mcmcfit, pars = c(thetaLabel, y0Label, vLabel, logpostLabel), separate_chains = T))

  #print(stan_diag(mcmcfit, information = "sample"))
  #print(stan_par(mcmcfit, par = c(thetaLabel, y0Label, vLabel, logpostLabel)))
  #print(stan_ess(mcmcfit, pars = c(thetaLabel, y0Label, vLabel, logpostLabel)))
  #print(stan_mcse(mcmcfit, pars = c(thetaLabel, y0Label, vLabel, logpostLabel)))
}


finalsolutionplot = function(mcmcList = NULL, xdata, ydata, meanSol = F)
{
  nruns = length(mcmcList)
  
  for(i in seq_len(nruns))
  {
    
    if(i == 1)
      plot(xdata, ydata, ylim=c(0,0.8), xlim=c(0, max(xdata)+5), main="Final Solution", xlab = "Time (Months)", ylab = "Predicted BCR-ABL level")
    
    #Set parameters
    mcmcData = as.matrix(mcmcList[[i]])
    colnames(mcmcData) = NULL
    thetamcmc = mcmcData[,8:15]
    y0mcmc = mcmcData[,16:20]
    vmcmc = mcmcData[,1]
    zmcmc = mcmcData[,2]
    loglikmcmc = mcmcData[,21]
    
    
    if(meanSol == T)
    {
      mcmcLen = dim(mcmcData)[1]
      
      times = seq(0, max(xdata) + 5, by = 0.1)
      R = matrix(nrow = mcmcLen, ncol = length(times))
      for(j in seq_len(mcmcLen))
      {
        theta = thetamcmc[j,]
        y0 = y0mcmc[j,]
        
        parameters = c(ax = theta[1], bx = theta[2], cx = theta[3],
                       dx = theta[4], ex = theta[5], ay = theta[6],
                       by = theta[7], ey = theta[8])
        
        state = c(x0 = y0[1], x1 = y0[2], x2 = y0[3], y0 = y0[4], y1 = y0[5])
        
        out = ode(y = state, times = times, func = desystem, parms = parameters)
        
        R[j,] = out[,6]/(out[,6] + 2*out[,4])
      }
      meanR = apply(R,2,mean)
      lines(times, meanR, col=i+1)
    }
    
    else
    {
      max_index = which.max(loglikmcmc)
      
      theta = thetamcmc[max_index,]
      y0 = y0mcmc[max_index,]
      
      print(theta)
      print(y0)
      print(loglikmcmc[max_index])
      
      parameters = c(ax = theta[1], bx = theta[2], cx = theta[3],
                     dx = theta[4], ex = theta[5], ay = theta[6],
                     by = theta[7], ey = theta[8])
      
      state = c(x0 = y0[1], x1 = y0[2], x2 = y0[3], y0 = y0[4], y1 = y0[5])
      
      times = seq(0, max(xdata) + 5, by = 0.1)
      out = ode(y = state, times = times, func = desystem, parms = parameters)
      
      R = out[,6]/(out[,6] + 2*out[,4])
      lines(out[,1], R, col=i+1)
    }
  }
}

PlotSampleSolutions = function(mcmcList = NULL, xdata, ydata, numSamples = 50, plotIndSamples = F)
{
  nruns = length(mcmcList)
  mcmcData = NULL
  for(i in seq_len(nruns))
    mcmcData = rbind(mcmcData, as.matrix(mcmcList[[i]]))
  
  colnames(mcmcData) = NULL
  thetamcmc = mcmcData[,8:15]
  y0mcmc = mcmcData[,16:20]
  vmcmc = mcmcData[,1]
  zmcmc = mcmcData[,2]
  loglikmcmc = mcmcData[,21]
  
  sampleInd = sample(seq_len(dim(mcmcData)[1]), numSamples)
  
  tempList = vector("list", numSamples)
  
  plot(xdata, ydata, ylim=c(0,1), xlim=c(0, max(xdata)+5), main="Sample Solutions", xlab = "Time (Months)", ylab = "Predicted BCR-ABL")
  
  times = seq(0, max(xdata) + 5, by = 0.1)
  
  for(i in seq_len(numSamples))
  {
    theta = thetamcmc[sampleInd[i],]
    y0 = y0mcmc[sampleInd[i],]
    
    parameters = c(ax = theta[1], bx = theta[2], cx = theta[3],
                   dx = theta[4], ex = theta[5], ay = theta[6],
                   by = theta[7], ey = theta[8])
    
    state = c(x0 = y0[1], x1 = y0[2], x2 = y0[3], y0 = y0[4], y1 = y0[5])
    
    out = ode(y = state, times = times, func = desystem, parms = parameters)
    
    tempList[[i]] = out
    
    R = out[,6]/(out[,6] + 2*out[,4])
    lines(times, R, col=8, lwd = 1)
  }
  points(xdata, ydata)
  
  if(plotIndSamples == T)
  {
    for(j in 1:5)
    {
      for(i in seq_len(numSamples))
      {
        tempout = tempList[[i]]
        if(i == 1)
          plot(tempout[,j+1], main = paste("Proportion of Variable",j), ylab = paste("Variable",j), xlab = "Time (Months)",  ylim = c(0,1), type = 'l', col=8, lwd=1)
        else
          lines(tempout[,j+1], col=8, lwd=1)
      }
    }
  }
  
  
}

IndividualSolutionPlot = function(mcmcList = NULL, xdata, ydata)
{
  nruns = length(mcmcList)
  
  for(j in 1:5)
  {
    for(i in seq_len(nruns))
    {
      #Set parameters
      mcmcData = as.matrix(mcmcList[[i]])
      colnames(mcmcData) = NULL
      thetamcmc = mcmcData[,8:15]
      y0mcmc = mcmcData[,16:20]
      vmcmc = mcmcData[,1]
      zmcmc = mcmcData[,2]
      loglikmcmc = mcmcData[,21]
      
      max_index = which.max(loglikmcmc)
      
      theta = thetamcmc[max_index,]
      y0 = y0mcmc[max_index,]
      
      parameters = c(ax = theta[1], bx = theta[2], cx = theta[3],
                     dx = theta[4], ex = theta[5], ay = theta[6],
                     by = theta[7], ey = theta[8])
      
      state = c(x0 = y0[1], x1 = y0[2], x2 = y0[3], y0 = y0[4], y1 = y0[5])
      
      times = seq(0, max(xdata) + 5, by = 0.1)
      out = ode(y = state, times = times, func = desystem, parms = parameters)

      if(i == 1)
        plot(times, out[,j+1], ylim=c(0,1), xlim=c(0, max(xdata)+5), ylab=paste("Proportion of Variable",j), main = paste("Trace of Variable",j), type='l')
      else
        lines(times, out[,j+1], col=i+1)
    }
  }
}


#plots the svd transform of the mcmc output in 3d
svdtransfromplot = function(mcmcList = NULL, mcmcData = NULL, angle = 30, dim = 3)
{
  if(is.null(mcmcData))
  {
    nruns = length(mcmcList)
    mcmcData = NULL
    for(i in seq_len(nruns))
      mcmcData = rbind(mcmcData, as.matrix(mcmcList[[i]]))
    
    colnames(mcmcData) = NULL
    thetamcmc = mcmcData[,8:15]
    y0mcmc = mcmcData[,16:20]
    vmcmc = mcmcData[,1]
    zmcmc = mcmcData[,2]
    loglikmcmc = mcmcData[,21]
  }
  else
  {
    mcmcData = as.matrix(mcmcData)
    colnames(mcmcData) = NULL
    thetamcmc = mcmcData[,8:15]
    y0mcmc = mcmcData[,16:20]
    vmcmc = mcmcData[,1]
    zmcmc = mcmcData[,2]
    loglikmcmc = mcmcData[,21]
  }
  
  m = cbind(thetamcmc,
            y0mcmc,
            vmcmc)
  
  n = dim(m)[1]
  k = dim(m)[2]
  
  numSamples = floor(n/10)
  
  decomp = svd(m)
  #sig = matrix(c(decomp$d[1], rep(0,n), decomp$d[2], rep(0, n), decomp$d[3], rep(0, n-3)), nrow=n)
  #transform = decomp$u %*% sig
  
  plot(decomp$d)
       
  #transform = m %*% decomp$v[,1:dim]
  transform = m %*% decomp$v[,c(1, k-1, k)]
  ind = sample(seq_len(n), numSamples)
  
  xrange = max(transform[ind,1]) - min(transform[ind,1])
  if(dim == 3)
    scatterplot3d(transform[ind,1],transform[ind,2],transform[ind,3], angle = 30, ylim = c(mean(transform[ind,2]) - xrange/2, mean(transform[ind,2]) + xrange/2), zlim = c(mean(transform[ind,3]) - xrange/2, mean(transform[ind,3]) + xrange/2))
  else
    plot(transform[ind,1], transform[ind,2], ylim = c(mean(transform[ind,2]) - xrange/2, mean(transform[ind,2]) + xrange/2))
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
    plot(PatientData[,2], PatientData[,3], ylim=c(0,1), xlab = "Time from beginning of treatment (Months)", ylab = "BCR-ABL level")
  
  return(PatientData)
}


dataPlot = function(fulldata, PatientInd = 25, overlay = F, type = 'l')
{
  tempdata = fulldata
  tempdata[,3] = as.numeric(tempdata[,3])
  tempdata = na.omit(tempdata)
  
  tempdata[,3] = tempdata[,3]/100
  count = 1
  for(PatientNumber in PatientInd)
  {
    ind = which(tempdata[,1] == PatientNumber)
    
    if(type != 'l')
      type = 'p'
    
    PatientData = tempdata[ind,]
    if(overlay == T)
    {
      if(count == 1)
        plot(PatientData[,2], PatientData[,3], ylim=c(0,1), xlim = c(0,60), xlab = "Time from beginning of treatment (Months)", ylab = "BCR-ABL level", type = type, col = count)
      else
        lines(PatientData[,2], PatientData[,3], col = count)
    }
    else
      plot(PatientData[,2], PatientData[,3], ylim=c(0,1), xlab = "Time from beginning of treatment (Months)", ylab = "BCR-ABL level", type = type)
    count = count + 1
  }  
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


priorSampling = function(n = 10, plotDensities = T, plotTraces = F, alpha = c(1,1,1,1,1))
{
  z = runif(n)
  #R0 = rbeta(n, 4,1)
  ytemp = rdirichlet(n, alpha)
  
  
  theta = matrix(runif(8*n), nrow = n)
  
  #ytemp = matrix(runif(5*n), nrow = n)
  #ysum = apply(ytemp, 1, sum)
  
  #ytemp = ytemp/ysum
  #y = ytemp
  #y[,3] = (ytemp[,3] + ytemp[,5])*(1-R0)/(1+R0)
  #y[,5] = (ytemp[,3] + ytemp[,5])*2*R0/(1+R0)
  
  y = z * ytemp
  
  R0 = y[,5]/(y[,5] + 2*y[,3])
  
  print(mean(R0))
  
  
  if(plotDensities == T)
  {
    par(mfrow = c(4,2))
    plot(density(theta[,1]), main="Prior density of theta 1")
    plot(density(theta[,2]), main="Prior density of theta 2")
    plot(density(theta[,3]), main="Prior density of theta 3")
    plot(density(theta[,4]), main="Prior density of theta 4")
    plot(density(theta[,5]), main="Prior density of theta 5")
    plot(density(theta[,6]), main="Prior density of theta 6")
    plot(density(theta[,7]), main="Prior density of theta 7")
    plot(density(theta[,8]), main="Prior density of theta 8")
    par(mfrow = c(3,2))
    plot(density(y[,1]), main="Prior density of IC for x0", xlab = "x0")
    plot(density(y[,2]), main="Prior density of IC for x1", xlab = "x1")
    plot(density(y[,3]), main="Prior density of IC for x2", xlab = "x2")
    plot(density(y[,4]), main="Prior density of IC for y0", xlab = "y0")
    plot(density(y[,5]), main="Prior density of IC for y1", xlab = "y1")
    plot(density(R0), main="Prior density of R0", xlab = "R0")
    par(mfrow = c(1,1))
    plot(density(z), main="Prior density of z")
  }
  
  if(plotTraces == T)
  {
    times = seq(0, 65, by = 0.1)
    R = matrix(nrow=n,ncol=length(times))
    tempList = vector("list", n)
    for(i in seq_len(n))
    {
      
      parameters = c(ax = theta[i,1], bx = theta[i,2], cx = theta[i,3],
                     dx = theta[i,4], ex = theta[i,5], ay = theta[i,6],
                     by = theta[i,7], ey = theta[i,8])
      
      state = c(x0 = y[i,1], x1 = y[i,2], x2 = y[i,3], y0 = y[i,4], y1 = y[i,5])
      
      out = ode(y = state, times = times, func = desystem, parms = parameters)
      tempList[[i]] = out
    }
    count = 0
    for(i in seq_len(n))
    {
      out = tempList[[i]]
      R = out[,6]/(out[,6] + 2*out[,4])
      if(R[2] < R[1])
        count = count + 1
      if(i == 1)
        plot(times, R, xlab="time", ylab="R", ylim = c(0,1), type='l', col = i, main = "Sample solutions with theta and ICs sampled from the priors")
      else
        lines(times, R, col = i)
    }
    print(count/n)
    for(i in seq_len(n))
    {
      out = tempList[[i]]
      if(i == 1)
        plot(times, out[,2], xlab="time", ylab="x0", ylim = c(0,1), col=i, type='l')
      else
        lines(times, out[,2], col = i)
    }
    for(i in seq_len(n))
    {
      out = tempList[[i]]
      if(i == 1)
        plot(times, out[,3], xlab="time", ylab="x1", ylim = c(0,1), type = 'l', col=i)
      else
        lines(times, out[,3], col = i)
    }
    for(i in seq_len(n))
    {
      out = tempList[[i]]
      if(i == 1)
        plot(times, out[,4], xlab="time", ylab="x2", ylim = c(0,1), type='l', col=i)
      else
        lines(times, out[,4], col = i)
    }
    for(i in seq_len(n))
    {
      out = tempList[[i]]
      if(i == 1)
        plot(times, out[,5], xlab="time", ylab="y0", ylim = c(0,1), type='l', col=i)
      else
        lines(times, out[,5], col = i)
    }
    for(i in seq_len(n))
    {
      out = tempList[[i]]
      if(i == 1)
        plot(times, out[,6], xlab="time", ylab="y1", ylim = c(0,1), type='l', col=i)
      else
        lines(times, out[,6], col = i)
    }
  }
}

#Generate Density and Trace plots from mcmc output
generatePlots = function(mcmcList = NULL, plotDensity = T, plotTrace = F, plotTheta = T, ploty0 = F, plotz = F, plotv = F, plotloglik = F, plotR0 = F)
{
  mcmcList = sampleMcmcList
  nruns = length(mcmcList)

  if(plotDensity == T)
  {
    if(plotTheta == T)
    {
      for(i in 1:8)
      {
        for(j in seq_len(nruns))
        {
          mcmcData = mcmcList[[j]]
          colnames(mcmcData) = NULL
          thetamcmc = mcmcData[,8:15]
          y0mcmc = mcmcData[,16:20]
          vmcmc = mcmcData[,1]
          zmcmc = mcmcData[,2]
          loglikmcmc = mcmcData[,21]
          
          if(j == 1)
            plot(density(thetamcmc[,i]), main=paste("Density of theta",i), col=j+1)
          else
            lines(density(thetamcmc[,i]), col=j+1)
        }
      }
    }
    if(ploty0 == T)
    {
      for(i in 1:5)
      {
        for(j in seq_len(nruns))
        {
          mcmcData = mcmcList[[j]]
          colnames(mcmcData) = NULL
          thetamcmc = mcmcData[,8:15]
          y0mcmc = mcmcData[,16:20]
          vmcmc = mcmcData[,1]
          zmcmc = mcmcData[,2]
          loglikmcmc = mcmcData[,21]
          
          if(j == 1)
            plot(density(y0mcmc[,i]), main=paste("Density of y0 for ode",i), col=j+1)
          else
            lines(density(y0mcmc[,i]), col=j+1)
        }
      }
    }
    if(plotR0 == T)
    {
      for(j in seq_len(nruns))
      {
        mcmcData = mcmcList[[j]]
        colnames(mcmcData) = NULL
        thetamcmc = mcmcData[,8:15]
        y0mcmc = mcmcData[,16:20]
        vmcmc = mcmcData[,1]
        zmcmc = mcmcData[,2]
        loglikmcmc = mcmcData[,21]
        R0mcmc = y0mcmc[,5]/(y0mcmc[,5] + 2*y0mcmc[,3])
        
        if(j == 1)
          plot(density(R0mcmc), main="Density of R0", col=j+1)
        else
          lines(density(R0mcmc), col=j+1)
      }
    }
    if(plotz == T)
    {
      for(j in seq_len(nruns))
      {
        mcmcData = mcmcList[[j]]
        colnames(mcmcData) = NULL
        thetamcmc = mcmcData[,8:15]
        y0mcmc = mcmcData[,16:20]
        vmcmc = mcmcData[,1]
        zmcmc = mcmcData[,2]
        loglikmcmc = mcmcData[,21]
        
        if(j == 1)
          plot(density(zmcmc), main="Density of z", col=j+1)
        else
          lines(density(zmcmc), col=j+1)
      }
    }
    if(plotv == T)
    {
      for(j in seq_len(nruns))
      {
        mcmcData = mcmcList[[j]]
        colnames(mcmcData) = NULL
        thetamcmc = mcmcData[,8:15]
        y0mcmc = mcmcData[,16:20]
        vmcmc = mcmcData[,1]
        zmcmc = mcmcData[,2]
        loglikmcmc = mcmcData[,21]
        
        if(j == 1)
          plot(density(vmcmc), main="Density of v", col=j+1)
        else
          lines(density(vmcmc), col=j+1)
      }
    }
    if(plotloglik == T)
    {
      for(j in seq_len(nruns))
      {
        mcmcData = mcmcList[[j]]
        colnames(mcmcData) = NULL
        thetamcmc = mcmcData[,8:15]
        y0mcmc = mcmcData[,16:20]
        vmcmc = mcmcData[,1]
        zmcmc = mcmcData[,2]
        loglikmcmc = mcmcData[,21]
        
        if(j == 1)
          plot(density(loglikmcmc), main="Density of Posterior Log-Likelihood", col=j+1)
        else
          lines(density(loglikmcmc), col=j+1)
      }
    }
  }
  if(plotTrace == T)
  {
    if(plotTheta == T)
    {
      for(i in 1:8)
      {
        for(j in seq_len(nruns))
        {
          mcmcData = mcmcList[[j]]
          colnames(mcmcData) = NULL
          thetamcmc = mcmcData[,8:15]
          y0mcmc = mcmcData[,16:20]
          vmcmc = mcmcData[,1]
          zmcmc = mcmcData[,2]
          loglikmcmc = mcmcData[,21]
          
          if(j == 1)
            plot(thetamcmc[,i], main=paste("Trace of theta",i), type='l', col=j+1)
          else
            lines(thetamcmc[,i], col=j+1)
        }
      }
      
    }
    if(ploty0 == T)
    {
      for(i in 1:5)
      {
        for(j in seq_len(nruns))
        {
          mcmcData = mcmcList[[j]]
          colnames(mcmcData) = NULL
          thetamcmc = mcmcData[,8:15]
          y0mcmc = mcmcData[,16:20]
          vmcmc = mcmcData[,1]
          zmcmc = mcmcData[,2]
          loglikmcmc = mcmcData[,21]
          
          if(j == 1)
            plot(y0mcmc[,i],main=paste("Trace of y0 for ode",i), type='l', col=j+1)
          else
            lines(y0mcmc[,i], col=j+1)
        }
      }
      
    }
    if(plotz == T)
    {
      for(j in seq_len(nruns))
      {
        mcmcData = mcmcList[[j]]
        colnames(mcmcData) = NULL
        thetamcmc = mcmcData[,8:15]
        y0mcmc = mcmcData[,16:20]
        vmcmc = mcmcData[,1]
        zmcmc = mcmcData[,2]
        loglikmcmc = mcmcData[,21]
        
        if(j == 1)
          plot(zmcmc, main="Trace of z", type='l', col = j+1)
        else
          lines(zmcmc, col = j+1)
      }
    }
    if(plotv == T)
    {
      for(j in seq_len(nruns))
      {
        mcmcData = mcmcList[[j]]
        colnames(mcmcData) = NULL
        thetamcmc = mcmcData[,8:15]
        y0mcmc = mcmcData[,16:20]
        vmcmc = mcmcData[,1]
        zmcmc = mcmcData[,2]
        loglikmcmc = mcmcData[,21]
        
        if(j == 1)
          plot(vmcmc,main="Trace of v", type='l', col = j+1)
        else
          lines(vmcmc, col = j+1)
      }
    }
    if(plotloglik == T)
    {
      for(j in seq_len(nruns))
      {
        mcmcData = mcmcList[[j]]
        colnames(mcmcData) = NULL
        thetamcmc = mcmcData[,8:15]
        y0mcmc = mcmcData[,16:20]
        vmcmc = mcmcData[,1]
        zmcmc = mcmcData[,2]
        loglikmcmc = mcmcData[,21]
        
        if(j == 1)
          plot(loglikmcmc,main="Trace of Posterior Log Likelihood", type='l', col = j+1)
        else
          lines(loglikmcmc, col = j+1)
      }
    }
  }
}
