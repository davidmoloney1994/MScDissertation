setwd('/Users/David/Dropbox/Oxford/Dissertation/')

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
            iter = 100, chains = 1)

print(fit)
get_inits(fit)

sampleMcmcList = generateListData(fit, nruns = n)
finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
generatePlots(mcmcList = sampleMcmcList, plotDensity = T, plotTheta = T, ploty0 = T, plotz = T, plotv = T)

svdtransfromplot(stanfit = fit, Phi = 0, Theta = 100)



#Multiple runs


inputdata = setRemissionData(data, showPlot = T)
inputdata = setPatientData(data, 25, showPlot = T)

standata = list(T = length(inputdata[,1]), 
                R = inputdata$BCR.ABL1.ABL1....,
                t0 = 0,
                ts = inputdata$Month[-1])

n = 8

fit <- stan(file = 'MScDissertation/odemcmc.stan', data = standata, 
              iter = 20000, chains = n, warmup = 2000)


sampleMcmcList = generateListData(fit, nruns = n)

par(mfrow = c(3,2))
finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1...., meanSol = T)
generatePlots(mcmcList = sampleMcmcList, plotDensity = T, plotTheta = T, ploty0 = T, plotz = T, plotv = T)
generatePlots(mcmcList = sampleMcmcList, plotDensity = F, plotTrace = T, plotTheta = F, ploty0 = F, plotz = F, plotv = F, plotloglik = T)
IndividualSolutionPlot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
PlotSampleSolutions(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1...., numSamples = 50)

#Write and read to file

#write.csv(extract(fit,permuted=T), file='50000relapse.csv', row.names=F)
mcmcdat = read.csv('50000relapse.csv', header=T)
sampleMcmcList = generateListData(mcmcData = mcmcdat)
finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
generatePlots(mcmcList = sampleMcmcList, plotDensity = T, plotTheta = T, ploty0 = T, plotz = T, plotv = T)
svdtransfromplot(mcmcData = mcmcdat, Phi = 0, Theta = 30)








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


finalsolutionplot = function(mcmcList = NULL, xdata, ydata, meanSol = F)
{
  nruns = length(mcmcList)
  
  for(i in seq_len(nruns))
  {
    
    if(i == 1)
      plot(xdata, ydata, ylim=c(0,1), xlim=c(0, max(xdata)+5), main="Final Solution")
    
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

PlotSampleSolutions = function(mcmcList = NULL, xdata, ydata, numSamples = 50)
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
  
  plot(xdata, ydata, ylim=c(0,1), xlim=c(0, max(xdata)+5), main="Sample Solutions")
  
  times = seq(0, max(xdata) + 5, by = 0.1)
  
  for(i in sampleInd)
  {
    theta = thetamcmc[i,]
    y0 = y0mcmc[i,]
    
    parameters = c(ax = theta[1], bx = theta[2], cx = theta[3],
                   dx = theta[4], ex = theta[5], ay = theta[6],
                   by = theta[7], ey = theta[8])
    
    state = c(x0 = y0[1], x1 = y0[2], x2 = y0[3], y0 = y0[4], y1 = y0[5])
    
    out = ode(y = state, times = times, func = desystem, parms = parameters)
    
    R = out[,6]/(out[,6] + 2*out[,4])
    lines(times, R, col=8, lwd = 1)
  }
  points(xdata, ydata)
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
    thetamcmc = mcmcData[,8:15]
    y0mcmc = mcmcData[,16:20]
    vmcmc = mcmcData[,1]
    zmcmc = mcmcData[,2]
    loglikmcmc = mcmcData[,21]
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
generatePlots = function(mcmcList = NULL, plotDensity = T, plotTrace = F, plotTheta = T, ploty0 = F, plotz = F, plotv = F, plotloglik = F)
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
