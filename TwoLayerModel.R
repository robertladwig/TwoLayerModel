rm(list = ls())
library(deSolve)
library(gridExtra)
library(ggplot2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

if (file.exists('output.txt')) 
  #Delete file if it exists
  file.remove('output.txt')

TwoLayer <- function(t, y, parms){
  dTe <-  Q / Ve * Tin -
    Q / Ve * y[1] +
    ((vt(t) * At) / Ve) * (y[2] - y[1]) +
    Jsw(t) / (rho * cp * H) * (10^4/0.2388) +
    (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt((Acoeff * (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t))))))) * (1 - Rl)) / (rho * cp * H) -
    (eps * sigma * (y[1] + 273)^4) /(rho * cp * H) -
    (c1 * Uw(t) * (y[1] - Tair(t))) / (rho * cp * H) -
    (Uw(t) * ((4.596 * exp((17.27 * y[1]) / (237.3 + y[1]))) - (Acoeff * (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t))))))) / (rho * cp * H)
  
  dTh <-  (vt(t) * At) / (Vh) * (y[1] - y[2])
  
  qin <- Q / Ve * Tin
  qout <- - Q / Ve * y[1] 
  mix_e <- ((vt(t) * At) / Ve) * (y[2] - y[1]) 
  mix_h <- (vt(t) * At) / (Vh) * (y[1] - y[2])
  sw <- Jsw(t) / (rho * cp * H) * (10^4/0.2388) 
  lw <- (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt((Acoeff * (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t))))))) * (1 - Rl)) / (rho * cp * H)
  water_lw <- - (eps * sigma * (y[1] + 273)^4) /(rho * cp * H)
  conv <- - (c1 * Uw(t) * (y[1] - Tair(t))) / (rho * cp * H)
  evap <- - (Uw(t) * ((4.596 * exp((17.27 * y[1]) / (237.3 + y[1]))) - (Acoeff * (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t))))))) / (rho * cp * H)
  
  write.table(matrix(c(qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap), nrow=1), 'output.txt', append = TRUE,
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(list(c(dTe, dTh)))
}


bound <- matrix(c(seq(1,12,1),
                  169, 274, 414, 552, 651, 684, 642, 537, 397, 259, 160, 127,
                  8.3, 9., 13.5,13.9,21.8,24.7,29.4,26.6,24.9,15.,9.7,6.6,
                  2.8,3.3,4.9,4.,5.3,7.8,11.8,11.5,7.7,6.8,6.5,2.4,
                  11.6,11.7,16.4,15.6,16.6,16.7,12.7,11.7,14.,12.9,14.8,11.6), nrow = 12, byrow = FALSE)
Ve <- 175000
Vh <- 75000
At <- 11000
Ht <- 3
As <- 25000
Tin <- 10
Q <- 7500
Rl <- 0.03 # reflection coefficient (generally small, 0.03)
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <- 4.9 * 10^(-3) # Stefan-Boltzmann constant in [J (m2 d K4)-1]
eps <- 0.97 # emissivity of water
rho <- 1000 # density
cp <- 4000 # specific heat
c1 <- 0.47 # Bowen's coefficient

parameters <- c(Ve = 175000,
                Vh = 75000,
                At = 11000,
                Ht = 3,
                As = 25000,
                Tin = 10,
                Q = 7500,
                Rl = 0.03, # reflection coefficient (generally small, 0.03)
                Acoeff = 0.6, # coefficient between 0.5 - 0.7)
                sigma = 11.7 * 10^(-8), # Stefan-Boltzmann constant in [cal (cm2 d K4)-1]
                eps = 0.97, # emissivity of water
                rho = 1000, # density
                cp = 4000, # specific heat
                c1 = 0.47 # Bowen's coefficient
)

H <- Ve/As
Et <- 7.07 * 10^(-4)  * H^(1.1505)
vto <- Et/Ht * (86400/10000)
vt_mix <- 1. # high mixing coefficient for overturn
bound <- as.data.frame(cbind(bound, c(rep(vt_mix,3),rep(vto,6),rep(vt_mix,3))))

colnames(bound) <- c('Month','Jsw','Tair','Dew','Uw','vt')
bound$Uw <- 19.0 + 0.95 * (bound$Uw)^2

#boundary <- rbind(bound,bound,bound)
boundary <- bound
boundary$Month <- seq(1, nrow(boundary),1)

Jsw <- approxfun(x = boundary$Month, y = boundary$Jsw, method = "linear", rule = 2)
Tair <- approxfun(x = boundary$Month, y = boundary$Tair, method = "linear", rule = 2)
Dew <- approxfun(x = boundary$Month, y = boundary$Dew, method = "linear", rule = 2)
Uw <- approxfun(x = boundary$Month, y = boundary$Uw, method = "linear", rule = 2)
vt <- approxfun(x = boundary$Month, y = boundary$vt, method = "linear", rule = 2)


times <- seq(from = 0, to = nrow(boundary), by = 1/30)
yini <- c(5,5)

out <- ode(times = times, y = yini, func = TwoLayer, parms = parameters, method = 'rk4')

plot(out[,1], out[,2], col = 'red', main = paste('vt_mix =',vt_mix))
lines(out[,1], out[,3], col = 'blue')

result <- data.frame('Time' = seq(1, nrow(out), 1),
                     'WT_epi' = out[,2], 'WT_hyp' = out[,3])
g1 <- ggplot(result) +
  geom_line(aes(x=Time, y=WT_epi, col='Epilimnion')) +
  geom_line(aes(x=(Time), y=WT_hyp, col='Hypolimnion')) +
  labs(x = 'Simulated Time', y = 'WT in deg C')

output <- read.table('output.txt')
# qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap
output <- data.frame('qin'=output[,1],'qout'=output[,2],'mix_e'=output[,3],'mix_h'=output[,4],
                     'sw'=output[,5],'lw'=output[,6],'water_lw'=output[,7],'conv'=output[,8],
                     'evap'=output[,9])
output$balance <- apply(output,1, sum)

g2 <- ggplot(output) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = qin, col = 'Qin')) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = qout, col = 'Qout')) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = mix_e, col = 'Mix_e')) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = mix_h, col = 'Mix_h')) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = sw, col = 'SW')) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = lw, col = 'LW')) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = water_lw, col = 'Reflec')) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = conv, col = 'Conv')) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = evap, col = 'Evap')) +
  geom_line(aes(x = seq(from = 1, to = nrow(output), by = 1),y = balance, col = 'Sum'), linetype = "dashed") +
  labs(x = 'Simulated Time', y = 'Flux in deg C/d')

#pdf('Simulation.pdf')
grid.arrange(g1, g2, ncol =1)
#dev.off()

#pdf('Simulation.pdf')
g3 <- grid.arrange(g1, g2, ncol =1)
ggsave(file='visual_result.png', g3, dpi = 300,width = 200,height = 180, units = 'mm')
#dev.off()