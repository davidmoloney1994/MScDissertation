setwd('/Users/David/Dropbox/Oxford/Dissertation/')
getwd()


library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(deSolve)


data = read.csv("datatest.csv", header=T, stringsAsFactors=F)

data[,3] = as.numeric(data[,3])
dat = na.omit(data)

dat[,3] = dat[,3]/100

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

#Remission data
luk_dat_remission <- list(T = length(remissionDat[,1]), 
                R = remissionDat$BCR.ABL1.ABL1....,
                t0 = 0,
                ts = remissionDat$Month[-1]) #Remove first time value. It's not used in the ode solver

#Patient 25, relapse
patient = 25
ind = which(dat[,1] == patient)
plot(dat[ind,2], dat[ind,3], ylim=c(0,1))

relapseDat25 = dat[which(dat[,1] == 25),]

luk_dat_relapse25 <- list(T = length(relapseDat25[,1]), 
                          R = relapseDat25$BCR.ABL1.ABL1....,
                          t0 = 0,
                          ts = relapseDat25$Month[-1]) #Remove first time value. It's not used in the ode solver

y0remissioninit = c(0.289, 0.009, 0.006, 0.04, 0.066)
thetaremissioninit = c(0.7,0.7,0.6,0.5,0.1,0.1,0.1,0.9)

y0relapseinit = c(0.27,0.219,0.05,0.077,0.162)
thetarelapseinit = c(0.1,0.8,0.2,0.8,0.2,0.9,0.4,0.8)

#Full remission data, no initialization
fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_remission, 
            iter = 1000, chains = 1)
#Full remission data, initialization for initial values
fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_remission, iter = 1000, chains = 1,
            init = list(list(y0temp=y0remissioninit/sum(y0remissioninit), z=sum(y0remissioninit))))
#Full remission data, full initialization
fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_remission, iter = 100, chains = 1,
            init = list(list(y0temp=y0remissioninit/sum(y0remissioninit), z=sum(y0remissioninit), theta=thetaremissioninit)))
#Patient 25, no initialization
fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_relapse25, 
            iter = 1000, chains = 1)
#Patient 25, initial value initialization
fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_relapse25, iter = 1000, chains = 1,
            init = list(list(y0temp=y0relapseinit/sum(y0relapseinit), z=sum(y0relapseinit))))
#Patient 25, full initialization
fit <- stan(file = 'MScDissertation/odemcmc.stan', data = luk_dat_relapse25, iter = 1000, chains = 1,
            init = list(list(y0temp=y0relapseinit/sum(y0relapseinit), z=sum(y0relapseinit), theta=thetarelapseinit)))


print(fit)
get_inits(fit)

la = extract(fit, permuted = T) # return a list of arrays 
thetamcmc = la$theta
y0mcmc = la$y0

n = length(la$theta[,1])
n=14

m = cbind(la$theta,
la$y0,
la$v)

decomp = svd(m)
plot(decomp$d)

decomp$d[1:2]
sig = matrix(c(decomp$d[1], rep(0,n), decomp$d[2], rep(0, n-2)), nrow=n)

plot(decomp$u %*% sig)

max_index = which.max(la$lp__)
max(la$lp__)

theta = thetamcmc[max_index,]
y0 = y0mcmc[max_index,]

#theta = apply(thetamcmc, 2, mean)
#y0 = apply(y0mcmc, 2, mean)

parameters = c(ax = theta[1], bx = theta[2], cx = theta[3],
               dx = theta[4], ex = theta[5], ay = theta[6],
               by = theta[7], ey = theta[8])

state = c(x0 = y0[1], x1 = y0[2], x2 = y0[3], y0 = y0[4], y1 = y0[5])

#parameters = c(ax=0.2,bx=0.8,cx=0.3,dx=0.8,ex=0.2,ay=0.9,by=0.3,ey=0.8)

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

times = seq(0, 70, by = 0.1)
out = ode(y = state, times = times, func = desystem, parms = parameters)
head(out)

R = out[,6]/(out[,6] + 2*out[,4])
lines(out[,1], R, col='red')

#likelihood plot
plot(extract(fit,permuted=F)[,1,21], t='l')
dimnames(extract(fit,permuted=F))
