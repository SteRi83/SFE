# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()


# settings
s0    = 100    # current stock price
k     = 100    # strike price
tau   = 1      # time to maturity in years
n.t   = 88     # number of time steps (used for stock simulation and averaging)
r     = 0.03   # annual interest rate
sigma = 0.2    # volatility
n.mc  = 10000  # number MC samples


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
    
    # geometric and arithmetic average
    s.avg.geom.i = exp(sum(log(st.i))/length(st.i))
    s.avg.arit.i = sum(st.i)/length(st.i)
    s.avg.i      = c(s.avg.arit.i, s.avg.geom.i) 
    
    # discounted average strike option payoff
    v.asc.i = exp(-r*tau)*pmax(0,  st.i[length(st.i)] - s.avg.i) 
    v.asp.i = exp(-r*tau)*pmax(0, -st.i[length(st.i)] + s.avg.i)

    # discounted average price option payoff
    v.apc.i = exp(-r*tau)*pmax(0,  s.avg.i - k)
    v.app.i = exp(-r*tau)*pmax(0, -s.avg.i + k)
    
    v.i     = c(v.asc.i, v.asp.i, v.apc.i, v.app.i)
    
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


AsianOutTable = function(value, s0, k, tau, n.t, r, sigma, n.mc){
  
  options(digits=3)
  
  cat(" -------------------------------------------------", "\n",
      "Value of Asian options by Monte Carlo", "\n",
      "-------------------------------------------------", "\n",
      " s0    = ", sprintf("%.2f", s0), "\n",
      " k     = ", sprintf("%.2f", k), "\n",
      " tau   = ", sprintf("%.2f", tau), "\n",
      " n.t   = ", sprintf("%.0f", n.t), "\n",
      " r     = ", sprintf("%.2f", r), "\n",
      " sigma = ", sprintf("%.2f", sigma), "\n",
      " n.nm  = ", sprintf("%.0f", n.mc), "\n",
      "-------------------------------------------------", "\n",
      "avg. strike|  arithm. avg.         geom. avg." , "\n",
      "-----------|-------------------------------------" , "\n",
      "      call |", sprintf("%.3f", value[1,1]), "+/-", sprintf("%.3f", value[2,1]),
      "   ",sprintf("%.3f", value[1,2]), "+/-", sprintf("%.3f", value[2,2]), "\n",
      "      put  |", sprintf("%.3f", value[1,3]), "+/-", sprintf("%.3f", value[2,3]),
      "   ",sprintf("%.3f", value[1,4]), "+/-", sprintf("%.3f", value[2,4]), "\n",
      "-----------|-------------------------------------" , "\n",
      " avg. price|  arithm. avg.         geom. avg." , "\n",
      "-----------|-------------------------------------" , "\n",
      "      call |", sprintf("%.3f", value[1,5]), "+/-", sprintf("%.3f", value[2,5]),
      "   ",sprintf("%.3f", value[1,6]), "+/-", sprintf("%.3f", value[2,6]), "\n",
      "      put  |", sprintf("%.3f", value[1,7]), "+/-", sprintf("%.3f", value[2,7]),
      "   ",sprintf("%.3f", value[1,8]), "+/-", sprintf("%.3f", value[2,8]), "\n")  
}

AsianOutTable(MCAsian(s0, k, tau, n.t, r, sigma, n.mc), s0, k, tau, n.t, r, sigma, n.mc)


# plot(t, st.i, type ="l", lwd = 1)
# lines(t, s.avg.geom.i*rep(1,length(t)) , type ="l", lwd = 2, col = "red")
# lines(t, s.avg.arit.i*rep(1,length(t)) , type ="l", lwd = 2, col = "blue")
# lines(t, test , type ="l", lwd = 2, col = "green")
