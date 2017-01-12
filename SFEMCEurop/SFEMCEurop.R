# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()


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


EuropOutTable = function(price){
  
  options(digits=3)
  
  cat(" Price of European options by Monte Carlo", "\n",
      "-----------------------------------------", "\n",
      "      call |", sprintf("%.3f", price[1,1]), "+/-", sprintf("%.3f", price[1,2]), "\n",
      "      put  |", sprintf("%.3f", price[2,1]), "+/-", sprintf("%.3f", price[2,2]), "\n")  
  
  return(0)
}

# Start the clock!
ptm <- proc.time()

EuropOutTable(MCEurop(S0 = 120, K = 130, tau = 0.6, nt = 1, r = 0.03, sigma = 0.2, nMC = 300000))

# Stop the clock
proc.time() - ptm
