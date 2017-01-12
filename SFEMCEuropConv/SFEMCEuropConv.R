rm(list = ls(all = TRUE))
graphics.off()


# input
S0    = 40
K     = 35
tau   = 1
r     = 0.03
sigma = 0.2


MCStock = function(S0, tau, nt, r, sigma, nMC){
  
  dt = tau/nt
  t  = seq(0, tau, by = dt)
  
  St = S0*exp((r-sigma^2/2)*t)*exp(replicate(nMC, c(0,cumsum(sigma*sqrt(dt)*rnorm(nt)))))
  
  return(St)
}


MCEurop = function(S0, K, tau, nt, r, sigma, nMC){
  
  St     = MCStock(S0, tau, nt, r, sigma, nMC)
  
  # discounted payoff
  C      = exp(-r*tau)*matrix(pmax(0,St[nrow(St),]-K), 1, ncol(St))
  P      = exp(-r*tau)*matrix(pmax(0,K-St[nrow(St),]), 1, ncol(St))
  V      = rbind(C, P)
   
  # MC estimate
  V_mc   = rowSums(V)/ncol(St)
  
  # 95% confidence intervall for MC-estimate
  conf95 = 1.96*1/sqrt(ncol(St))*sqrt(1/(ncol(St)-1)*rowSums((V-V_mc)^2))
 
  result = cbind(V_mc, conf95)
  
  return(result)
  
}


BSCall = function(S, K, tau, r, sigma) {

  d1 = (log(S/K) + (r + (sigma^2/2)) * tau)/(sigma * sqrt(tau))
  d2 = d1 - (sigma * sqrt(tau))
  C  = S * pnorm(d1) - K * exp(-r * tau) * pnorm(d2)

  return(C)
}

k   = seq(2, 6, 0.5)
nMC = 10^k 

VCall = matrix(0,length(nMC),2)

for (i in 1:length(nMC)){
  
  Call = MCEurop(S0, K, tau, nt = 1, r, sigma, nMC = nMC[i])

  # MC estimate  
  VCall[i,1] = Call[1,1]
  
  # 95% confidence 
  VCall[i,2] = Call[1,2]
  
}

BSC = BSCall(S0, K, tau, r, sigma)


par(mar=c(5,5,2,1)) 
plot(nMC, VCall[,1], log = "x", type = "p", pch = 16, ylim = c(min(VCall[,1]-VCall[,2]), max(VCall[,1]+VCall[,2])),
     ylab = "Value of European Call", xlab = "Number of MC samples", cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
lines(nMC, matrix(BSC,length(nMC),1)[,1], type = "l", lwd = 2)
arrows(nMC, VCall[,1]-VCall[,2], nMC, VCall[,1]+VCall[,2], code=3, angle=90, length=0.1, lwd = 2)
legend(10^4, max(VCall[,1]+VCall[,2]), legend = c("MC-Estimation with 95% CI", "Black-Scholes-Model"), lty = c(0,1), pch = c(16,NA), cex = 1.2)

