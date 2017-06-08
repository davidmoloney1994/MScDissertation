data = read.csv("./Dropbox/Oxford/Dissertation/datatest.csv", header=T, stringsAsFactors=F)

library(deSolve)

data[,3] = as.numeric(data[,3])
dat = na.omit(data)

dat[,3] = dat[,3]/100

patient = 3
ind = which(dat[,1] == patient)

plot(dat[ind,2], dat[ind,3])

head(dat)
relapse = NULL

for(i in 1:(length(dat[,1]) - 1))
{
  if(dat[i+1,3] > dat[i,3] + 0.05 && dat[i,1] == dat[i+1,1])
    relapse = c(relapse, dat[i,1])
}

relapse = unique(relapse)
remission = setdiff(1:69,relapse)

remissionDat = dat[which(dat[,1] %in% remission),]
relapseDat = dat[which(dat[,1] %in% relapse),]

remissionDat = aggregate(BCR.ABL1.ABL1.... ~ Month, data = remissionDat, mean)

plot(remissionDat$Month, remissionDat$BCR.ABL1.ABL1....)

for(i in remission)
{
  patient = i
  ind = which(dat[,1] == patient)
  
  plot(dat[ind,2], dat[ind,3])
}

for(i in relapse)
{
  patient = i
  ind = which(dat[,1] == patient)
  
  plot(dat[ind,2], dat[ind,3])
}




theta = runif(8,0,1)
y0 = runif(5,0,0.2)
sigma = runif(0,1)


parameters = c(ax = 0.5, bx = 0.5, cx = 0.5,
               dx = 0.5, ex = 0.5, ay = 0.5,
               by = 0.5, ey = 0.5)

state = c(x0 = 0.1, x1 = 0.1, x2 = 0.1, y0 = 0.1, y1 = 0.1)

parameters = c(ax = theta[1], bx = theta[2], cx = theta[3],
               dx = theta[4], ex = theta[5], ay = theta[6],
               by = theta[7], ey = theta[8])

state = c(x0 = y0[1], x1 = y0[2], x2 = y0[3], y0 = y0[4], y1 = y0[5])


desystem <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx0 <- ax * x0 * (1 - x0 - x1 - x2 - y0 - y1) - bx * x0
    dx1 <- bx * x0 + cx * x1 * (1 - x0 - x1 - x2 - y0 - y1) - dx * x1
    dx2 <- dx * x1 - ex * x2
    dy0 <- ay * y0 * (1 - x0 - x1 - x2 - y0 - y1) - by * y0
    dy1 <- by * y0 - ey * y1
    list(c(dx0, dx1, dx2, dy0, dy1))
  })
}

times <- seq(0, 70, by = 1)




out <- ode(y = state, times = times, func = desystem, parms = parameters)
head(out)

R = out[,6]/(out[,6] + 2*out[,4])

lines(out[,1], R, col='red')




