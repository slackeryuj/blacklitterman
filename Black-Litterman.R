
# install.packages("fPortfolio")
# install.packages("BLCOP")
# install.packages("mnormt")
# install.packages("tseries")
# install.packages("zoo")
# install.packages("PerformanceAnalytics")
# install.packages("Rsolnp")
# install.packages("quadprog")
# install.packages("rgl")
# install.packages("rugarch")
# install.packages("fGarch")
# install.packages("VineCopula")
# install.packages("BLModel")

library(fPortfolio)
library(BLCOP)
library(mnormt)
library(BLModel)


library(tseries)
library(zoo)
library(PerformanceAnalytics)
library(Rsolnp)
library(quadprog)
library(rgl)
library(rugarch)
library(VineCopula)
library(MASS)
library(fGarch)
library(rmgarch)  #install.packages("rmgarch")
library(ggplot2)
library(GGally)


rm(list=ls())

setwd("C:/Users/Jiayang/Desktop/OR 664/")


# United States (GDP: 20.49 trillion) SPY
# China (GDP: 13.4 trillion) GXC
# Japan: (GDP: 4.97 trillion) EWJ
# Germany: (GDP: 4.00 trillion) EWG
# France: (GDP: 2.78 trillion) EWQ

# Brazil: (GDP: 1.87 trillion) EWZ
# Canada: (GDP: 1.71 trillion) EWC

# United States
# China
# Japan
# Germany
# India
# United Kingdom
# France
# Italy
# Brazil
# Canada


# United Kingdom: (GDP: 2.83 trillion) EWU
# Italy: (GDP: 2.07 trillion) EWI
# India: (GDP: 2.72 trillion) PIN

## specify the tikcers universe
spComp <- read.delim("temp_ticker.txt", header = FALSE) 


## specify time period
dateStart <- "2007-03-30"               
dateEnd <- "2020-5-28"

symbols <- spComp[, 1]
nAss <- length(symbols)


z <- get.hist.quote(instrument = as.character(symbols[1]), start = dateStart,
                    end = dateEnd,
                    quote = "Adjusted",
                    provider = c("yahoo"), 
                    retclass = c("zoo", "ts"), quiet = FALSE, drop = FALSE)
dimnames(z)[[2]] <- as.character(symbols[1])


if(length(symbols)>1){
  for(i in 2:nAss) {
    ## display progress by showing the current iteration step
    cat("Downloading ", as.character(symbols[i]), " out of ", nAss , "\n")
    
    result <- try(x <- get.hist.quote(instrument = as.character(symbols[i]), start = dateStart,
                                      end = dateEnd,
                                      quote = "Adjusted",
                                      provider = c("yahoo"), 
                                      retclass = c("zoo", "ts"), quiet = FALSE, drop = FALSE))
    if(class(result) == "try-error") {
      next
    }
    else {
      dimnames(x)[[2]] <- as.character(symbols[i])
      
      ## merge with already downloaded data to get assets on same dates 
      z <- merge(z, x)                      
      
    }
  }
}

BL.poseterior <- function(pickMatrix, qv,confidences, mu, tau = 0.5, sigma, kappa = 0) {
  
  P <- as.matrix(pickMatrix)
  sigma <- as.matrix(sigma)
  mu <- matrix(mu,ncol = 1)
  qv <- matrix(qv,ncol = 1)
  if (kappa == 0) {
    if (length(confidences) > 1) 
      omega <- diag(1/confidences)
    else omega <- matrix(1/confidences, 1, 1)
  } else {
    omega <- kappa * tcrossprod(P %*% sigma, P)
    omegaInv <- solve(omega)
  }
  sigmaInv <- solve(sigma)
  temp <- tcrossprod(sigma, P)
  postMu <- mu + tau * temp %*% solve(tau * P %*% temp + omega, qv - P %*% mu)
  postterior_result = list()
  postterior_result[[1]] <- as.numeric(postMu)
  names(postterior_result[[1]] ) <- colnames(country_return)
  postSigma <- (1 + tau) * sigma - tau^2 * temp %*% solve(tau *P %*% temp + omega, P %*% sigma)
  postterior_result[[2]] <- postSigma
  names(postterior_result)=c("PosteriorMean", "PosteriorVariance")
  return(postterior_result)
}



data.monthly <- z[ endpoints(z, on="months", k=1), ]

country_return <- diff(log(data.monthly), lag=1)

country_return=country_return[,-c(5,7,8)]

country_return = country_return[index(country_return)<='2019-12-31',]

# country_return <- na.omit(country_return)

par(mfrow=c(3,3))

acfPlot(country_return)

ggpairs(country_return)


for (i in 1:7){
  qqnorm(country_return[,i], main=colnames(country_return)[i])
  qqline(country_return[,i])
}



nAsset = length(colnames(country_return))

view_matrix <- read.csv("view_matrix.csv", header = T, row.names = 1)
priorMeans <- matrix(NA,1,nAsset)

yy = 1

final_weight = matrix(NA,48, nAsset)
portfolio_ret = matrix(NA,48, 4)
colnames(portfolio_ret) = c("Optimal", "EW", "Optimal-Prior", "Optimal-Prior-Simulation")

final_weight_prior = matrix(NA,48, nAsset)
final_weight_prior_simu = matrix(NA,48, nAsset)

for (i in 1: 48){
  
  country_return_use = country_return[i:(i+104),]
  
  regression_result = CAPMList(country_return_use, marketIndex = country_return_use$SPY, riskFree = 0, regFunc = "rlm")
  
  for (j in 1:nAsset) {
    priorMeans[,j] = mean(regression_result$alphas[j]+regression_result$betas[j]*country_return_use$SPY)
  }
  
  # specify i.i.d. model for the univariate time series
  ugarch_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
                            variance.model = list(model = "sGARCH", garchOrder = c(1,1)))
  
  # specify DCC model
  dcc_spec <- dccspec(uspec = multispec(replicate(ugarch_spec, n = nAsset)),
                      VAR = TRUE, lag = nAsset,
                      model = "DCC", dccOrder = c(1,1))
  garchdcc_fit <- dccfit(dcc_spec, data = country_return_use, solver = "nlminb")
  # garchdcc_fit
  
  dcc.fcst= dccforecast(garchdcc_fit, n.ahead=1)
  
  priorSigma <- rcov(dcc.fcst)
  
  pickMatrix <- matrix(0, 2, nAsset)
  pickMatrix[1,which(view_matrix[,yy]>0)] = -view_matrix[which(view_matrix[,yy]>0)[1],yy] #overvalued assets - underperform S&P 500
  pickMatrix[1,1] = 1
  pickMatrix[2,which(view_matrix[,yy]<0)] = -view_matrix[which(view_matrix[,yy]<0)[1],yy] #undervalued assets - outperform S&P 500
  pickMatrix[2,1] = -1
  
  # pickMatrix[1,which(view_matrix[,yy]>0)] = view_matrix[which(view_matrix[,yy]>0)[1],yy] #undervalued assets - outperform S&P 500
  # pickMatrix[1,1] = -1
  # pickMatrix[2,which(view_matrix[,yy]<0)] = view_matrix[which(view_matrix[,yy]<0)[1],yy] #overvalued assets - underperform S&P 500
  # pickMatrix[2,1] = 1

  qv =  c(0.1/12*mean(regression_result$betas[which(view_matrix[,yy]>0)]), 0.1/12*mean(regression_result$betas[which(view_matrix[,yy]<0)]))
  
  
  country_return_posterior = BL.poseterior(pickMatrix, qv,confidences=NA, priorMeans, tau = 1/105, priorSigma[[1]][,,1], kappa = 1)
  
  set.seed(123)
  tsData = timeSeries(mvrnorm(n = 1000, mu = country_return_posterior$PosteriorMean, Sigma = country_return_posterior$PosteriorVariance))
  tsData_prior_simu = timeSeries(mvrnorm(n = 1000, mu = priorMeans, Sigma = cov(country_return_use)))
  
  tsData_prior = timeSeries(country_return_use)
  
  
  
  frontierSpec <- portfolioSpec()
  
  
  # setSolver(frontierSpec) <-"solveRglpk.CVAR"
  
  # setAlpha( frontierSpec ) <-0.01
  setNFrontierPoints( frontierSpec ) <-50
  setRiskFreeRate(frontierSpec)<- 0.025
  
  # solveRglpk.CVAR(tsData, frontierSpec, constraints ="LongOnly")
  
  # Constraint = c("minW[1  :nAssets]=0", "maxW[1:nAssets]=0.3")
  
  frontier1g <- portfolioFrontier( data =tsData , spec = frontierSpec , constraints = "LongOnly")
  muTab1g <- getTargetReturn( frontier1g@portfolio )
  weightTab1g <- getWeights( frontier1g@portfolio )
  

  ret_tgt <- mean(tsData)
  
  
  ind <-1
  while ( muTab1g[ind ,1] < min(ret_tgt)){
    ind <-ind +1
  } 
  final_weight[i,] = weightTab1g[ind ,]
  
  prior_optimalPort <- portfolioFrontier( data =tsData_prior , spec = frontierSpec , constraints = "LongOnly")
  mu_prior <- getTargetReturn( prior_optimalPort@portfolio )
  weight_prior <- getWeights( prior_optimalPort@portfolio )
  ind2 <-1
  while ( mu_prior[ind2 ,1] < min(ret_tgt)){
    ind2 <-ind2 +1
  } 
  final_weight_prior[i,] = weight_prior[ind ,]
  
  prior_optimalPort_simu <- portfolioFrontier( data =tsData_prior_simu , spec = frontierSpec , constraints = "LongOnly")
  mu_prior_simu <- getTargetReturn( prior_optimalPort_simu@portfolio )
  weight_prior_simu <- getWeights( prior_optimalPort_simu@portfolio )
  ind3 <-1
  while ( mu_prior_simu[ind3 ,1] < min(ret_tgt)){
    ind3 <-ind3 +1
  } 
  final_weight_prior_simu[i,] = weight_prior_simu[ind ,]
  
  
  
  
  # frontier1g <- tangencyPortfolio( data =tsData , spec = frontierSpec , constraints = "LongOnly")
  # 
  # 
  # final_weight[i,] = frontier1g@portfolio@portfolio$weights
  # 
  # prior_optimalPort <- tangencyPortfolio( data =tsData_prior , spec = frontierSpec , constraints = "LongOnly")
  # 
  # 
  # final_weight_prior[i,] = prior_optimalPort@portfolio@portfolio$weights
  
  portfolio_ret[i,1] =  final_weight[i,]%*%t(country_return[i+105,])
  portfolio_ret[i,2] =  rep(1/nAsset,nAsset)%*%t(country_return[i+105,])
  portfolio_ret[i,3] =  final_weight_prior[i,]%*%t(country_return[i+105,])
  portfolio_ret[i,4] =  final_weight_prior_simu[i,]%*%t(country_return[i+105,])
  
  if (i%%12==0) {yy = yy+1}
  
  cat("Finish", as.character(i), " out of ", 48 , "\n")
}

out_ret=country_return[106:(length(country_return[,1])),]

Port_Value =  zoo(ts(100*cumprod(portfolio_ret[,1]+1)), index(out_ret))
BMK_Value =   zoo(ts(100*cumprod(portfolio_ret[,2]+1)), index(out_ret))
Port_Value_Prior =  zoo(ts(100*cumprod(portfolio_ret[,3]+1)), index(out_ret)) 
Port_Value_Prior_Simu =  zoo(ts(100*cumprod(portfolio_ret[,4]+1)), index(out_ret)) 

plot(Port_Value, col = 'blue', ylim=c(0,200), main='Out Of Sample Result vs. Benchmark', xlab='Date', ylab = 'Portfolio Values ($100 Deposit)')
lines(Port_Value_Prior, col = 'red')
lines(Port_Value_Prior_Simu, col = 'green')
legend(index(Port_Value)[1],200,c("Proppsed Model","Benchmark - Historical", "Benchmark - Historical Simulation"),lty=c(1,1,1),col=c("blue","red", "green"), cex = 0.8)

((coredata(Port_Value[48])/100)^(1/4)-1)/(StdDev(coredata(portfolio_ret[,1]))*sqrt(12))

((coredata(Port_Value_Prior[48])/100)^(1/4)-1)/(StdDev(coredata(portfolio_ret[,3]))*sqrt(12))

((coredata(Port_Value_Prior_Simu[48])/100)^(1/4)-1)/(StdDev(coredata(portfolio_ret[,4]))*sqrt(12))

