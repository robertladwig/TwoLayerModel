rm(list = ls())
library(deSolve)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(bvpSolve)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

if (file.exists('output.txt')) 
  #Delete file if it exists
  file.remove('output.txt')

OneLayer <- function(t, y, parms){
  with(as.list(c(y, parms)),{
  
  eair <- (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t)))) # air vapor pressure
  esat <- 4.596 * exp((17.27 * Tair(t)) / (237.3 + Tair(t))) # saturation vapor pressure
  RH <- eair/esat *100 # relative humidity
  # epilimnion water temperature change per time unit
  dT <-  Q / V * Tin -              # inflow heat
    Q / V * y +                   # outflow heat
    + As/(V * rho * cp) * (
      Jsw(t)  + # shortwave radiation
        (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
        (eps * sigma * (y + 273)^4)  - # backscattering longwave radiation from the lake
        (c1 * Uw(t) * (y - Tair(t))) - # convection
        (Uw(t) * ((esat) - (eair))) )# evaporation
  
  qin <- Q * rho * cp * Tin /As
  qout <- - Q * rho * cp * y /As
  sw <- Jsw(t) 
  lw <- (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) 
  water_lw <- - (eps * sigma * (y + 273)^4)
  conv <- - (c1 * Uw(t) * (y - Tair(t))) 
  evap <- - (Uw(t) * ((esat) - (eair)))
  Rh <- RH
  
  write.table(matrix(c(qin, qout, sw, lw, water_lw, conv, evap, Rh, dT, y, Jsw(t), Tair(t), Dew(t), Uw(t)), nrow=1), 'output.txt', append = TRUE,
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(list(c(dT)))
  })
}


bound <- matrix(c(seq(1,12,1),
                  169, 274, 414, 552, 651, 684, 642, 537, 397, 259, 160, 127,
                  8.3, 9., 13.5,13.9,21.8,24.7,29.4,26.6,24.9,15.,9.7,6.6,
                  2.8,3.3,4.9,4.,5.3,7.8,11.8,11.5,7.7,6.8,6.5,2.4,
                  11.6,11.7,16.4,15.6,16.6,16.7,12.7,11.7,14.,12.9,14.8,11.6), nrow = 12, byrow = FALSE)
V <- 250000 * 10^6
As <- 25000 * 10^4
Tin <- 10
Q <- 7500 * 10^6
Rl <- 0.03 # reflection coefficient (generally small, 0.03)
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <- 11.7 * 10^(-8)#4.9 * 10^(-3) # Stefan-Boltzmann constant in [J (m2 d K4)-1]
eps <- 0.97 # emissivity of water
rho <- 0.9982 # density
cp <- 0.99 # specific heat
c1 <- 0.47 # Bowen's coefficient

parameters <- list(V = V,As = As,Tin =Tin,Q = Q,Rl = Rl,Acoeff = Acoeff,sigma = sigma,
                   eps = eps,rho = rho,cp = cp,c1 =c1)

colnames(bound) <- c('Month','Jsw','Tair','Dew','Uw')
bound <- as.data.frame(bound)
bound$Uw <- 19.0 + 0.95 * (bound$Uw * 1000/3600)^2

#boundary <- rbind(bound,bound,bound)
boundary <- bound
boundary$Month <- seq(1, nrow(boundary),1)
#boundary$Jsw <- boundary$Jsw * 100

Jsw <- approxfun(x = boundary$Month, y = boundary$Jsw, method = "linear", rule = 2)
Tair <- approxfun(x = boundary$Month, y = boundary$Tair, method = "linear", rule = 2)
Dew <- approxfun(x = boundary$Month, y = boundary$Dew, method = "linear", rule = 2)
Uw <- approxfun(x = boundary$Month, y = boundary$Uw, method = "linear", rule = 2)


times <- seq(from = 1, to = nrow(boundary), by = 1/30)
yini <- c(y = 8)

out <- ode(times = times, y = yini, func = OneLayer, parms = parameters, method = 'rk4')
diagnostics(out)

plot(out[,1], out[,2], col = 'red')

result <- data.frame('Time' = seq(1, nrow(out), 1),
                     'WT' = out[,2])
g1 <- ggplot(result) +
  geom_line(aes(x=Time, y=WT, col='1D Model')) +
  labs(x = 'Simulated Time', y = 'WT in deg C')  +
  theme_bw()

output <- read.table('output.txt')
sim.output <- data.frame('dTe' = output[,9],'Te'=output[,10], 'Jsw' = output[,11], 'Tair' = output[,12],
                         'Dew' = output[,13], 'Uw' = output[,14])
# qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap
diag.output <- data.frame('qin'=output[,1],'qout'=output[,2],
                     'sw'=output[,3],'lw'=output[,4],'water_lw'=output[,5],'conv'=output[,6],
                     'evap'=output[,7], 'Rh' = output[,8])
diag.output$balance <- apply(diag.output[,-c(8)],1, sum)

g2 <- ggplot(diag.output) +
  geom_line(aes(x = seq(from = 1, to = nrow(diag.output), by = 1),y = qin, col = 'Qin')) +
  geom_line(aes(x = seq(from = 1, to = nrow(diag.output), by = 1),y = qout, col = 'Qout')) +
  geom_line(aes(x = seq(from = 1, to = nrow(diag.output), by = 1),y = sw, col = 'SW')) +
  geom_line(aes(x = seq(from = 1, to = nrow(diag.output), by = 1),y = lw, col = 'LW')) +
  geom_line(aes(x = seq(from = 1, to = nrow(diag.output), by = 1),y = water_lw, col = 'Reflec')) +
  geom_line(aes(x = seq(from = 1, to = nrow(diag.output), by = 1),y = conv, col = 'Conv')) +
  geom_line(aes(x = seq(from = 1, to = nrow(diag.output), by = 1),y = evap, col = 'Evap')) +
  geom_line(aes(x = seq(from = 1, to = nrow(diag.output), by = 1),y = balance, col = 'Sum'), linetype = "dashed") +
  scale_colour_brewer("Colors in Set1", palette="Set3") +
  labs(x = 'Simulated Time', y = 'Flux in cal/(cm2 d)')  +
  theme_bw()

#pdf('Simulation.pdf')
grid.arrange(g1, g2, ncol =1)
#dev.off()

#pdf('Simulation.pdf')
g3 <- grid.arrange(g1, g2, ncol =1)
ggsave(file='1L_visual_result.png', g3, dpi = 300,width = 200,height = 180, units = 'mm')
#dev.off()