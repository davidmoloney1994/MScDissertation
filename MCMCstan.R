setwd('/Users/David/Dropbox/Oxford/Dissertation/')

#Various pakages to be installed and options to be set for stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(deSolve)
library(plot3D)
library(rBeta2009)
library(scatterplot3d)
library(shinystan)
library(parallel)
library(MASS)
library(maxmatching)
library(igraph)
library(pdist)
library(bmk)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(plotrix)

#Read in the data
data = read.csv("datatest.csv", header=T, stringsAsFactors=F)

#Set the data to be used
inputdata = setRemissionData(data, showPlot = T)
inputdata = setPatientData(data, 25, showPlot = F)

#set data for stan to use
standata = list(T = length(inputdata[,1]), 
                R = inputdata$BCR.ABL1.ABL1....,
                t0 = 0,
                ts = inputdata$Month[-1],
                alpha = c(0.9,0.9,0.4,0.9,1.5))

#number of mcmc chains
n = 4

#Runs mcmc and saves in 'fit'
fit <- stan(file = 'MScDissertation/odemcmc.stan', data = standata, 
              iter = 30000, chains = n, warmup = 3000)

#Save the fit to file
#saveRDS(fit, file = "10Remission30000Normalfit.RData")


#Create list of mcmc output from 'fit' for use in functions
sampleMcmcList = generateListData(samplefit, nruns = 4)

par(mfrow = c(1,2))
finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1...., meanSol = T)
generatePlots(mcmcList = sampleMcmcList, plotDensity = T, plotTheta = T, ploty0 = T, plotz = F, plotv = T, plotR0 = T, plotloglik = T)
generatePlots(mcmcList = sampleMcmcList, plotDensity = F, plotTrace = T, plotTheta = F, ploty0 = F, plotz = F, plotv = F, plotloglik = T)
IndividualSolutionPlot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
PlotSampleSolutions(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1...., numSamples = 30, plotIndSamples = F)
svdtransfromplot(mcmcList  = sampleMcmcList, angle = 20, dim=2)
generatePlots(mcmcList = sampleMcmcList, plotDensity = F, plotTrace = T, plotTheta = T, ploty0 = T, plotz = F, plotv = T, plotloglik = T)

samplefit = readRDS("25Relapse30000Normalfit.RData")
samplefit = readRDS("10Remission30000Normalfit.RData")
generatediagnostics(mcmcfit = samplefit, plotTrace = T, plotACF = T, plotDensity = F, plotTheta = F, ploty0 = F, plotv = F, plotloglik = T)
print(samplefit)

######################################################
#Prior Elicitation
######################################################
priorSampling(n = 50, plotTraces = T, plotDensities = F, alpha = c(0.9,0.9,0.4,0.9,1.5))
priorSampling(n = 1000000, plotTraces = F, plotDensities = T, alpha = c(0.9,0.9,0.4,0.9,1.5)) 


v = rgamma(1000000,0.25,5)
plot(density(v), main = "", xlab = "v")

median(sqrtv)

plot(density(sqrt(v)), xlim = c(0,1), main="density of sampled standard deviation", xlab = "sqrt(v)")

######################################################
#EDA
######################################################
remmissionInd = c(10, 23)
relapseInd = c(25, 47)

par(mfrow = c(1,2))
#plot data for patients in PatientInd
dataPlot(data, PatientInd = remmissionInd, overlay = F, type = 'l')

######################################################
#LDA
######################################################
PatientIDs = 1:69
PatientIDs = PatientIDs[-66]
relapse_patient_IDs = c(7,25,47,26,1,68)

ldafit = PatientLDA(PatientIDs = PatientIDs, theta_ind = 1:8, y0_ind = NULL, v_ind = NULL, num_samples = 10000, relapse_IDs = relapse_patient_IDs)
#saveRDS(ldafit, file = "ldafit.RData")

ldafit = readRDS("ldafit.RData")

xlim = c(-3,2)
ylim = c(0,22)
px = sort(ldafit$scaling)
py = rep(0,8)
lx.buf = 0.2
lx = seq(xlim[1]+lx.buf,xlim[2]-lx.buf,len=length(px))
ly = 5

par(xaxs='i',yaxs='i',mar=c(5,1,1,1))
plot(NA,xlim=xlim,ylim=ylim,axes=F,ann=F)
axis(1)

segments(px,py,lx,ly, lty = 2)
points(px,py,pch=16,xpd=NA, col=2)
text(lx,ly, c(expression(w[7]) , expression(w[1]), expression(w[8]), expression(w[3]), expression(w[6]), expression(w[5]), expression(w[4]), expression(w[2])),pos=3)
lab = order(as.numeric(ldafit$scaling))
expression(paste(theta[lab[1]]))

dev.off()

######################################################
#Clustering
######################################################
PatientIDs = c(10,25)
PatientIDs = PatientIDs[-66]

DistanceMatrix = DifferenceMatrix(PatientNumbers = PatientIDs, paired_dist = F, avg_dist = F, theta_ind = 1:8) 
#saveRDS(DistanceMatrix, file = "ThetaALL_DistMatrix.RData")
CILengthMatrix = DifferenceMatrix(PatientNumbers = PatientIDs, paired_dist = F, avg_dist = T, theta_ind = 1:8, ci = T)
max(CILengthMatrix)


DistanceMatrix = readRDS("EuclideanMean_DistMatrix.RData")
image(DistanceMatrix)
hc = hclust(as.dist(DistanceMatrix), "complete")
dend = as.dendrogram(hc, hang=0.05)
dend1 = color_branches(dend, k = 4, col = c(4,2,3,1))
plot(dend1, ylim = c(0.3,1.6))
abline(h=1.3, lty= 2)


DistanceMatrix = readRDS("Theta7_DistMatrix.RData")
mds = cmdscale(as.dist(DistanceMatrix), eig = T, k=2)
#mds = sammon(as.dist(DistanceMatrix), magic=0.51)
vec = rep(0.7,68)
vec[c(7, 25, 47, 26, 1, 67)] = 1
cols = rep(1,68)
cols[c(7, 25, 47, 26, 1, 67)] = cols[c(7, 25, 47, 26, 1, 67)] + 1
plot(mds$points[,1], mds$points[,2], xlab="Coordinate 1", ylab="Coordinate 2",	type="n")
text(mds$points[,1], mds$points[,2], labels = row.names(mds$points), cex=vec, col=cols)

par(mfrow = c(1,1))
par(mar=c(5.1, 4.1, 0.5, 0.5))
dev.off()
PatientDiagnostics(PatientID = 29, plotTrace = F, plotACF = F, plotDensity = F, plotTheta = F, ploty0 = F, plotv = F, plotloglik = F, plotsol = T)
legend(30,0.65, c("data point","predicted BCR-ABL"), pch = c(1,NA), lty = c(NA, 1), col = c(1,2))







######################################################
##################   Functions   #####################
######################################################

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

#Takes a stan fit object and returns a list of the data
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

clusterGroupAnalysis = function(hclustObject, no_groups = 4, selected_group = 1)
{
  groups = cutree(hclustObject, no_groups)
  group = PatientIDs[groups == selected_group]
  for(i in group)
    PatientDiagnostics(PatientID = i, plotsol = T)
  
  print("Patents:")
  print(group)
}

#Generates various diagnostic plots for stan fit object 'mcmcfit'. It can plot a trace/density of theta/y0/v/logpost
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
    print(stan_trace(mcmcfit, pars = c(thetaLabel, y0Label, vLabel, logpostLabel), inc_warmup = TRUE, nrow = 2) + xlim(0, 10000)) 
  
  if(plotACF == T)
    print(stan_ac(mcmcfit, pars = c(thetaLabel, y0Label, vLabel, logpostLabel))+ ylab("Autocorrelation")) 
  
  if(plotDensity == T)
    print(stan_dens(mcmcfit, pars = c(thetaLabel, y0Label, vLabel, logpostLabel), separate_chains = F, nrow=4,ncol=2))
}

#Calls the function 'generatediagnostics' for patient 'PatientID
#If plotsol = T then the predicted BCR-ABL is plotted also
PatientDiagnostics = function(PatientID = 1, plotTrace = F, plotACF = F, plotDensity = F, plotTheta = F, ploty0 = F, plotv = F, plotloglik = F, plotsol = F)
{
  samplefit = readRDS(paste("PatientData/",PatientID, "_10000fit.RData", sep=""))
  
  if(plotsol == T)
  {
    data = read.csv("datatest.csv", header=T, stringsAsFactors=F)
    inputdata = setPatientData(data, PatientID, showPlot = F)
    sampleMcmcList = generateListData(samplefit, nruns = 4)
    finalsolutionplot(mcmcList = sampleMcmcList, xdata = inputdata$Month, ydata = inputdata$BCR.ABL1.ABL1....)
  }
  generatediagnostics(mcmcfit = samplefit, plotTrace = plotTrace, plotACF = plotACF, plotDensity = plotDensity, plotTheta = plotTheta, ploty0 = ploty0, plotv = plotv, plotloglik = plotloglik)
}

#Calculates the total KS distance summed over theta parameter 'theta_ind'
Total_CDF_Distance = function(Patients = c(1,2), n = 4, theta_ind = 1:8)
{
  samplefit1 = readRDS(paste("PatientData/",Patients[1], "_10000fit.RData", sep=""))
  samplefit2 = readRDS(paste("PatientData/",Patients[2], "_10000fit.RData", sep=""))

  data1 = extract(samplefit1)
  data2 = extract(samplefit2)
  
  mat1 = data1$theta
  mat2 = data2$theta
  
  totalDiff = 0
  ind = seq(0,1,0.001)
  for(k in theta_ind)
  {
    cdf1 = ecdf(mat1[,k])
    cdf2 = ecdf(mat2[,k])
    totalDiff = totalDiff + max(abs(cdf1(ind) - cdf2(ind)))
  }
  return(totalDiff)
}

#Calculates the approximate mean euclidea distance between samples from patients 'Patients'
#If CI is TRUE, it instead calculates the width of the confidence interval
Average_Distance = function(Patients = c(1,2), num_samples = 100000, type = "mean", CI = F)
{
  samplefit1 = readRDS(paste("PatientData/",Patients[1], "_10000fit.RData", sep=""))
  samplefit2 = readRDS(paste("PatientData/",Patients[2], "_10000fit.RData", sep=""))
  
  data1 = extract(samplefit1)
  data2 = extract(samplefit2)
  
  mat1 = data1$theta
  mat2 = data2$theta
  
  ind1 = sample(seq_len(dim(mat1)[1]), num_samples, replace = T)
  ind2 = sample(seq_len(dim(mat2)[1]), num_samples, replace = T)
  
  eucl = numeric(num_samples)
  for(i in seq_len(num_samples))
    eucl[i] = dist(rbind(mat1[ind1[i], ], mat2[ind2[i], ]))

  if(CI == T)
  {
    sig = sd(eucl)
    tsig = qt(1 - 0.05/2, df = num_samples - 1)
    return(2 * tsig * sig / sqrt(num_samples))
  }
  
  else
  {
    if(type == "mean")
      return(mean(eucl))
    
    if(type == "max")
      return(max(eucl))
  }
}

best.col.perm<-function(y,yp)
{
  
  #returns a matrix yp[,v] derived from yp by permuting its cols
  #the permutation v is chosen to minimise sum(|y-yp[,v]|)
  
  n=dim(y)[2]; if (dim(yp)[2]!=n) stop('different number of cols');  
  p=dim(y)[1]; if (dim(yp)[1]!=p) stop('different number of rows');
  
  #The matrix dm[] measures the closeness of match between columns
  #cols of y are vertices 1:n, cols of yp are n+1 to 2n
  
  dm=matrix(0,2*n,2*n)
  
  #A 0 in this matrix will mean no edge. There is no edge between
  #vertices corresponding to cols of the same matrix.
  #The closeness of the match is a num between 1 and 2
  #2=perfect match, cols are equal, 1=cols are perfect opposites.
  
  db=as.matrix(pdist(X=t(y),Y=t(yp)))
  a=1.01*max(db)-db #maximum matching so we need similarities not distances
  b=matrix(0,n,n)
  dm=rbind(cbind(b,a),cbind(t(a),b))
  
  G=graph_from_adjacency_matrix(dm, mode="undirected", weighted=TRUE)
  G=set_vertex_attr(G, name="type", value=c(rep(TRUE,n),rep(FALSE,n)))
  
  mt=maxmatching(G, weighted = TRUE, maxcardinality = TRUE)
  v=mt$matching[(1:n)]-n
  return(yp[,v])
  
}

mdist = function(y,yp)
{
  n=dim(y)[1];
  if (dim(yp)[1]!=n) stop('different number of rows');
  
  p=dim(y)[2];
  if (dim(yp)[2]!=p) stop('different number of columns');
  
  ypp=t(best.col.perm(t(y),t(yp)))
  
  d=sum(sqrt(apply((y-ypp)^2,1,sum)))
  return(d)
}

Paired_Distance = function(Patients = c(1,2), num_samples = 1000)
{
  samplefit1 = readRDS(paste("PatientData/",Patients[1], "_10000fit.RData", sep=""))
  samplefit2 = readRDS(paste("PatientData/",Patients[2], "_10000fit.RData", sep=""))
  
  data1 = extract(samplefit1)
  data2 = extract(samplefit2)
  
  mat1 = data1$theta
  mat2 = data2$theta
  
  ind1 = sample(seq_len(dim(mat1)[1]), num_samples, replace = F)
  ind2 = sample(seq_len(dim(mat2)[1]), num_samples, replace = F)
  
  return(mdist(mat1[ind1,], mat2[ind2,]))
}

#Calculates the distance matrix for one of the 3 defined distances
#avg_dist = T   => use mean distance between samples
#avg_dist = F   => use KS distance with theta parameter 'theta_ind'
DifferenceMatrix = function(PatientNumbers = 1:69, n = 4, theta_ind = 1:8, avg_dist = F, paired_dist = F, ci = F)
{
  n = length(PatientNumbers)
  grid = combn(PatientNumbers, 2)
  
  no_cores = detectCores() - 1
  cl = makeCluster(no_cores)
  clusterExport(cl, "extract")
  clusterExport(cl, "mdist")
  clusterExport(cl, "best.col.perm")
  clusterExport(cl, "pdist")
  clusterExport(cl, "graph_from_adjacency_matrix")
  clusterExport(cl, "set_vertex_attr")
  clusterExport(cl, "maxmatching")
  
  if(avg_dist == T && paired_dist == T)
    stop("Only one distance required")
  
  if(ci == T && avg_dist == F)
    stop("Confidence interval length only for avg distance")
  {
    if(avg_dist == T)
      temp = parApply(cl, grid, 2, Average_Distance, num_samples = 100000, type = "mean", CI = ci)
    
    else if(paired_dist == T)
      temp = parApply(cl, grid, 2, Paired_Distance, num_samples = 1000)
    
    else
      temp = parApply(cl, grid, 2, Total_CDF_Distance, n = 4, theta_ind = theta_ind)
  }
  
  stopCluster(cl)
  
  mat = matrix(0, nrow = n, ncol = n)
  mat[lower.tri(mat)] = temp
  mat = mat + t(mat)
  colnames(mat) = PatientNumbers
  return(mat)
}

#Returns LDA fit object with 'num_samples' samples taken from each patient in 'PatientIDs'
#Patients 'relapse_IDs' are categorised as relapse
PatientLDA = function(PatientIDs = 1:69, theta_ind = 1:8, y0_ind = 1:5, v_ind = 1, num_samples = 1000, relapse_IDs = c(7,25,47))
{
  group = rep(as.numeric(PatientIDs %in% relapse_IDs) + 1, each = num_samples)
  samples = SampleFromPost(PatientIDs = PatientIDs, theta_ind = theta_ind, y0_ind = y0_ind, v_ind = v_ind, num_samples = num_samples)
  
  ldafit = lda(x = samples, grouping = group)
  return(ldafit)
}

#Plots the mode predicted BCR-ABL. If ci=T then confience intervals are plotted
#meanSol = T   =>   the mean BCR-ABL at each time is plotted.  !!!! Long time to run !!!!
#indSol = T  =>  the mode solution for each chain is plotted
#input is a stan ouput list as returned by 'generateListData'
#xdata and ydata also need to be set and inputed as returned by 'setPatientData'
finalsolutionplot = function(mcmcList = NULL, xdata, ydata, meanSol = F, indSol = F, ci = F)
{
  nruns = length(mcmcList)
  mcmcData = NULL
  if(indSol == T)
  {
    for(i in seq_len(nruns))
    {
      if(i == 1)
        plot(xdata, ydata, ylim=c(0,0.8), xlim=c(0, max(xdata)+5), xlab = "Time (Months)", ylab = "Predicted BCR-ABL level")
      
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
  else
  {
    for(i in seq_len(nruns))
      mcmcData = rbind(mcmcData, as.matrix(mcmcList[[i]]))
    
    colnames(mcmcData) = NULL
    thetamcmc = mcmcData[,8:15]
    y0mcmc = mcmcData[,16:20]
    vmcmc = mcmcData[,1]
    zmcmc = mcmcData[,2]
    loglikmcmc = mcmcData[,21]
    
    max_index = which.max(loglikmcmc)
    
    theta = thetamcmc[max_index,]
    y0 = y0mcmc[max_index,]
    v = vmcmc[max_index]
    sig = sqrt(v)
    
    
    parameters = c(ax = theta[1], bx = theta[2], cx = theta[3],
                   dx = theta[4], ex = theta[5], ay = theta[6],
                   by = theta[7], ey = theta[8])
    
    state = c(x0 = y0[1], x1 = y0[2], x2 = y0[3], y0 = y0[4], y1 = y0[5])
    
    times = seq(0, max(xdata) + 5, by = 0.1)
    out = ode(y = state, times = times, func = desystem, parms = parameters)

    ind = match(xdata, times)
    
    R = out[,6]/(out[,6] + 2*out[,4])
    
    L = R[ind] - 1.96 * sig
    U = R[ind] + 1.96 * sig
    
    if(ci == T)
    {
      plotCI(xdata, R[ind], ui=U, li=L, xlim=c(0, max(xdata)+5), ylim=c(min(0,L), max(U)+0.05) , xlab = "Time (Months)", ylab = "Predicted BCR-ABL level", pch=NA)
      points(xdata, ydata)
      
      lines(out[,1], R, col=2)
      abline(h=0)
      legend(25,0.65, c("data point","predicted BCR-ABL"), pch = c(1,NA), lty = c(NA, 1), col = c(1,2))
    }
    else
    {
      plot(xdata, ydata, xlim=c(0, max(xdata)+5), ylim=c(0, max(ydata)+0.05) , xlab = "Time (Months)", ylab = "Predicted BCR-ABL level")
      lines(out[,1], R, col=2)
    }
  }
}

#Generates samples from the posterior for a given patient and return them in a matrix.
SampleFromPost = function(PatientIDs = 1:69, theta_ind = 1:8, y0_ind = 1:5, v_ind = 1, num_samples = 1000)
{
  sampleMat = NULL
  for(i in PatientIDs)
  {
    samplefit = readRDS(paste("PatientData/",i, "_10000fit.RData", sep=""))
    data = extract(samplefit)
    theta = data$theta
    y0 = data$y0
    v = data$v
    
    samp_ind = sample(dim(theta)[1], num_samples)
    
    if(is.null(v_ind))
      PatientSample = cbind(theta[samp_ind, theta_ind], y0[samp_ind, y0_ind])
    else
      PatientSample = cbind(theta[samp_ind, thetaInd], y0[samp_ind, y0_ind], v[samp_ind])
      
    sampleMat = rbind(sampleMat, PatientSample)
  }
  return(as.matrix(sampleMat))
}

#Function to plot BCR-ABL level plots of posterior samples. Similar input to 'finalsolutionplot'
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
  
  plot(xdata, ydata, ylim=c(0,0.8), xlim=c(0, max(xdata)+5), xlab = "Time (Months)", ylab = "Predicted BCR-ABL")
  
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

#Plots the individual solutions y1,...y5 to the system of ODEs
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

#Plots the data for given patients 'PatientInd'
#If overlay = T, then all patient in PatientInd are plotted in one plot
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
      plot(PatientData[,2], PatientData[,3], ylim=c(0,1), xlab = "Time (Months)", ylab = "BCR-ABL level", type = type)
    count = count + 1
  }  
}

#Sets aggrigated remission data
#Not used in the report
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

#Samples from the prior.
#plotDensities = T   =>   Density plots of the priors are returned
#plotTrace = T   =>   prior BCR-ABL solution plots are returned
#alpha is the coefficient for the unscaled Dirichlet in the prior
priorSampling = function(n = 10, plotDensities = T, plotTraces = F, alpha = c(1,1,1,1,1))
{
  z = runif(n)
  ytemp = rdirichlet(n, alpha)
  
  theta = matrix(runif(8*n), nrow = n)
  
  y = z * ytemp
  
  R0 = y[,5]/(y[,5] + 2*y[,3])
  
  print(mean(R0))
  finalVals = numeric(n)
  
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
    par(mfrow = c(1,1))
    plot(density(R0), xlab = "R0")
    plot(density(z), main="Prior density of z")
    par(mfrow = c(3,2))
    plot(density(y[,1]), main = "", xlab = "y1")
    plot(density(y[,2]), main = "", xlab = "y2")
    plot(density(y[,3]), main = "", xlab = "y3")
    plot(density(y[,4]), main = "", xlab = "y4")
    plot(density(y[,5]), main = "", xlab = "y5")
  }
  
  if(plotTraces == T)
  {
    times = seq(0, 65, by = 0.1)
    L = length(times)
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
      finalVals[i] = R[L]
      if(R[2] < R[1])
        count = count + 1
      if(i == 1)
        plot(times, R, xlab="Time (Months)", ylab="BCR-ABL level", ylim = c(0,1), type='l', col = i)
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
    
    print("Percent less than 0.05:")
    print(sum(finalVals < 0.05)*100/n)
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
