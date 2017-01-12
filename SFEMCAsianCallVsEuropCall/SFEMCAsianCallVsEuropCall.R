# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()


MCPath = function(s0, tau, n.t, r, sigma){
  
  dt = tau/n.t
  t  = seq(0, tau, by = dt)
  
  st = s0*exp((r-sigma^2/2)*t + sigma*sqrt(dt)*c(0, cumsum(rnorm(n.t))))
  
  return(st)

}

MCAsian = function(s0, k, tau, n.t, r, sigma, n.mc){
  
  v.sum         = 0
  v.squared.sum = 0
  
  for (i in 1:n.mc){
    
    st.i = MCPath(s0, tau, n.t, r, sigma)
    
    # arithmetic average
    s.avg.i = sum(st.i)/length(st.i)

    # discounted average strike option payoff
    v.asc.i = exp(-r*tau)*pmax(0,  st.i[length(st.i)] - s.avg.i) 

    # discounted average price option payoff
    v.apc.i = exp(-r*tau)*pmax(0,  s.avg.i - k)

    v.i     = c(v.asc.i, v.apc.i)
    
    # calculate sums for MC estimate and confidence intervall
    v.sum         = v.sum + v.i
    v.squared.sum = v.squared.sum + v.i^2
  }  

  # MC estimate
  v.mc = v.sum/n.mc
  
  # 95% confidence intervall
  conf95 = 1.96*1/sqrt(n.mc)*sqrt(1/(n.mc-1)*(v.squared.sum - 2*v.mc*v.sum + n.mc*v.mc^2))

  return(rbind(v.mc, conf95))
}


BSEuropCall = function(S, K, tau, r, sigma) {
  
  d1 = (log(S/K) + (r + (sigma^2/2)) * tau)/(sigma * sqrt(tau))
  d2 = d1 - (sigma * sqrt(tau))
  C  = S * pnorm(d1) - K * exp(-r * tau) * pnorm(d2)
  
  return(C)
}

s0        = seq(25, 45, by = 0.5)
AsianCall = rbind(rep(0, length(s0)),rep(0, length(s0)))
EuropCall = rep(0, length(s0))

for(i in 1:length(s0)){
  
  k     = 35
  tau   = 1
  n.t   = 88
  r     = 0.03
  sigma = 0.2
  n.mc  = 10^5
  print(i)
  
  AsianCall[1:2,i] = MCAsian(s0[i], k, tau, n.t, r, sigma, n.mc)[1:2,2]
  EuropCall[i] = BSEuropCall(s0[i], k, tau, r, sigma)
  
}
pdf("asianvseurop.pdf", width=8, height=5)
par(mar=c(5,5,2,1)) 
plot (s0, AsianCall[1,], type ="l", lwd = 2, col = "black", ylim = c(min(AsianCall, EuropCall), max(AsianCall, EuropCall)),
      xlab = "Asset spot price",
      ylab = "Option value", cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
lines(s0, EuropCall, type ="l", lwd = 2, lty = 5,  col = "black")
legend(min(s0), max(EuropCall), legend = c("Black-Scholes European call", "MC Asian average price call"), lwd = 2, lty = c(5,1), cex = 1.2)
dev.off()
