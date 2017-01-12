# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# s0  - current stock price
# tau - time to maturity in years
# n.t - number of time steps


# function for generating geometric brownian motion
MCPath = function(s0, tau, n.t, r, sigma){
  
  dt = tau/n.t
  t  = seq(0, tau, by = dt)
  
  st = s0*exp((r-sigma^2/2)*t + sigma*sqrt(dt)*c(0, cumsum(rnorm(n.t))))
  
  return(st)

}

MCAsian = function(s0, k, tau, n.t, r, sigma){
  
  exit   = FALSE
  n.loop = 0
  
  while (exit == FALSE && n.loop < 1e4){
    
    n.loop  = n.loop + 1

    st.i    = MCPath(s0, tau, n.t, r, sigma)
    s.avg.i = sum(st.i)/length(st.i)

    if (st.i[length(st.i)] > 1.05*s.avg.i && s.avg.i > 1.05*k){   # factor to prevent overlaying text 
      PayoffPlot(st.i, k, tau, n.t, s.avg.i)
      exit = TRUE
    } 
  }     
  

}


PayoffPlot = function(st.i, k, tau, n.t, s.avg.arit.i){
  
  dt = tau/n.t
  t  = seq(0, tau, by = dt)
  
  f.wind = 1.7
  f.line = 1.15
  f.text = 1.2
  f.arro = 1.1 
  
  s.T = st.i[length(st.i)]
  
  # path plot
  par(mar=c(5,5,2,1)) 
  plot(t, st.i,
       type ="l", 
       lwd  = 2, 
       xlim = c(0, f.wind*tau),
       ylim = c(min(k,st.i), max(k, st.i)),
       xaxs ="i", 
       cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2,
       main = "Asian Call payoff", 
       xlab = "Time t", 
       ylab = expression("Stockprice "*"S"["t"]))
  # S_T 
  lines(c(tau, f.line*tau ), s.T*c(1,1) , lty = 5, lwd = 2, col = "black")
  text(f.text*tau, s.T, expression("S"["T"]), col = "black", cex = 1.2)
  # S_avg
  lines(c(0, f.line*tau), s.avg.arit.i*c(1,1) , lty = 5, lwd = 2, col = "black")
  text(f.text*tau, s.avg.arit.i, expression("  S"["avg"]), col = "black", cex = 1.2)
  # K
  lines(c(0, f.line*tau), k*c(1,1)            , lty = 5, lwd = 2, col = "black")
  text(f.text*tau, k, "K", col = "black", cex = 1.2)
  # S_avg <--> K avg. price payoff
  arrows(f.arro*tau, s.T, f.arro*tau, s.avg.arit.i, code = 3, length=0.1)
  text(1.4*tau, 0.5*(k+s.avg.arit.i), expression("(S"["avg"]*"  -  "*"K)"*"   avg. price payoff"), cex = 1.2)
  # S_t <--> s_avg avg. strike payoff
  arrows(f.arro*tau, k, f.arro*tau, s.avg.arit.i, code = 3, length=0.1)
  text(1.4*tau, 0.5*(s.T+s.avg.arit.i), expression("(S"["T"]*"  -  "*"S"["avg"]*")"*"   avg. strike payoff"), cex = 1.2)
}


MCAsian(50, 55, 1, 500, 0.03, 0.2)


