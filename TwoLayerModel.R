######
library(deSolve)
rigidode <- function(t,y,parms){
  dy1 <- -2 * y[2] * y[3]
  dy2 <- 1.25 * y[1] * y[3]
  dy3 <- -0.5 * y[1] * y[2]
  list(c(dy1, dy2, dy3))
}

yini <- c(y1=1,y2=0,y3=0.9)
times <- seq(from = 0, to = 20, by = 0.01)
out <- ode(times = times, y = yini, func = rigidode, parms = NULL, method = 'rk4')
######

######
library(bvpSolve)
fun <- function(t, y, pars) { 
  dy1 <- y[2]
  dy2 <- - 3 * p * y[1] / (p+t*t)^2
  return(list(c(dy1,
                dy2))) }
p    <- 1e-5

init <- c(-0.1 / sqrt(p+0.01), NA)
end  <- c( 1 / sqrt(p+0.01), NA)

# Solve bvp
sol  <- bvpcol(yini = init, yend = end, 
               x = seq(-0.1, 0.1, by = 0.001), func = fun)
plot(sol, which = 1)
######

######
library(deSolve)
working <- structure(list(datetime = structure(c(1185915600, 1185919200, 
                                                 1185922800, 1185926400, 1185930000, 1185933600, 1185937200, 1185940800, 
                                                 1185944400, 1185948000, 1185951600), class = c("POSIXct", "POSIXt"
                                                 ), tzone = "UTC"), p = c(0, 0, 0, 1.1, 0.5, 0.7, 0, 0, 1.3, 0, 
                                                                          0), e = c(0.15, 0.14, 0.13, 0.21, 0.15, 0.1, 0.049, 0, 0, 0, 
                                                                                    0), qsim = c(-1.44604436552566, NA, NA, NA, NA, NA, NA, NA, NA, 
                                                                                                 NA, NA)), .Names = c("datetime", "p", "e", "qsim"), row.names = c(NA, 11L), 
                     class = "data.frame")
gradfun <- function(t,y,parms) {
  pcp <- working$p[pmax(1,ceiling(t))]
  pet <- working$e[pmax(1,ceiling(t))]
  list(y^2*((pcp-pet)/y),NULL)
}

gradfun(0,working$qsim[1],1)   ## test
ode1 <- ode(y=c(qsim=working$qsim[1]),func=gradfun,
            time=seq(0,nrow(working),length.out=101),
            parms=NULL,method="rk4")
plot(ode1)
######

######
library(deSolve)

exponential=function(t,state,parameters){ with(as.list( c(state,parameters)), {
  
  #Aux. Var.
  fX2 = pmax(0,1-(1-(d2/r12)*(X2/K2)))
  fX1 = X1/(X1+k1); 
  
  # equations (ODE)
  dX1 = C-((d1)*(X1))-(r12)*(X2)*fX2*fX1 # differential equaion
  dX2 = r12*(X2)*fX2*fX1-((d2)*(X2))
  
  return(list(c(dX1, dX2, K2)))
})
}

# -- RUN INFORMATION

# Set Initial Values and Simulation Time
state = c(X1=2,X2=0.01,K2= 10)
times=0:100

# Assign Parameter Values
parameters = c(d1=0.001, d2=0.008, r12=0.3,C=0.5,k1= 0.001)

for (i in 1:length(times)){
  out= ode(y=state,times=times,func=exponential,parms=parameters)
}
######

######
library(deSolve)
eqRG = function(tm, state, parms)
{
  with(as.list(c(tm, state, parms)),
       {
         a1 = parms[["a1"]]
         a2 = parms[["a2"]]
         a3 = parms[["a3"]]
         b1 = parms[["b1"]]
         b2 = parms[["b2"]]
         b3 = parms[["b3"]]
         c1 = parms[["c1"]]
         c2 = parms[["c2"]]
         c3 = parms[["c3"]]
         dy1 = a1(tm) * y2 - b1(tm) * y3 - c1(tm) * y1
         dy2 = a2(tm) * y1 - b2(tm) * y1 - c2(tm) * y3
         dy3 = a3(tm) * y2 - b3(tm) * y3 - c3(tm) * y1
         return(list(c(dy1, dy2, dy3)))
       })
}
tm = seq(0, 10, len = 100)
state = c(y1 = 1, y2 = 0.5, y3 = 0.02)
a1 = b1 = c1 = approxfun(tm, - tm/10)
a2 = b2 = c2 = approxfun(tm, tm * 2)
a3 = b3 = c3 = approxfun(tm, sin(tm / 20))
P = list(a1 = a1, a2 = a2, a3 = a3,
         b1 = a1, b2 = b2, b3 = b3,
         c1 = c1, c2 = c2, c3 = c3)
sol = ode(y = state, times = tm, parms = P, func = eqRG)
plot(sol)
######

######
## Forcing function data
Flux <- matrix(ncol=2,byrow=TRUE,data=c(
  1, 0.654, 11, 0.167,   21, 0.060, 41, 0.070, 73,0.277, 83,0.186,
  93,0.140,103, 0.255,  113, 0.231,123, 0.309,133,1.127,143,1.923,
  153,1.091,163,1.001,  173, 1.691,183, 1.404,194,1.226,204,0.767,
  214, 0.893,224,0.737, 234,0.772,244, 0.726,254,0.624,264,0.439,
  274,0.168,284 ,0.280, 294,0.202,304, 0.193,315,0.286,325,0.599,
  335, 1.889,345, 0.996,355,0.681,365,1.135))

parms <- c(k=0.01)

times <- 1:365

## the model
sediment <- function( t, O2, k) 
  list (c(Depo(t) - k * O2), depo = Depo(t))

# the forcing functions; rule = 2 avoids NaNs in interpolation
Depo <- approxfun(x = Flux[,1], y = Flux[,2], method = "linear", rule = 2)

Out <- ode(times = times, func = sediment, y = c(O2 = 63), parms = parms)
######

rm(list = ls())
library(deSolve)

TwoLayer <- function(t, y, parms){
  dTe <-  Q / Ve * Tin -
    Q / Ve * y[1] +
    ((vt(t) * At) / Ve) * (y[2] - y[1]) +
    Jsw(t) / (rho * cp * Ht) +
    (rho * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt((Acoeff * (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t))))))) * (1 - Rl)) / (rho * cp * Ht) -
    (eps * sigma * (y[1] + 273)^4) /(rho * cp * Ht) -
    (c1 * Uw(t) * (y[1] - Tair(t))) / (rho * cp * Ht) -
    (Uw(t) * ((4.596 * exp((17.27 * y[1]) / (237.3 + y[1]))) - (Acoeff * (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t))))))) / (rho * cp * Ht)
  dTh <-  (vt(t) * At) / (Vh) * (y[1] - y[2])
  
  return(list(c(dTe, dTh)))
}

boundary <- matrix(c(seq(1,12,1),
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
H <- Ve/As
Et <- 7.07 * 10^(-4)  * H^(1.1505)
vt <- Et/Ht * (86400/10000)
vt_mix <- 1 # high mixing coefficient for overturn
Rl <- 0.03 # reflection coefficient (generally small, 0.03)
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <- 11.7 * 10^(-8) # Stefan-Boltzmann constant in [cal (cm2 d K4)-1]
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

boundary <- as.data.frame(cbind(boundary, c(rep(vt_mix,4),rep(vt,5),rep(vt_mix,3))))
colnames(boundary) <- c('Month','Jsw','Tair','Dew','Uw','vt')
boundary$Uw <- 19.0 + 0.95 * (boundary$Uw)^2

Jsw <- approxfun(x = boundary$Month, y = boundary$Jsw, method = "linear", rule = 2)
Tair <- approxfun(x = boundary$Month, y = boundary$Tair, method = "linear", rule = 2)
Dew <- approxfun(x = boundary$Month, y = boundary$Dew, method = "linear", rule = 2)
Uw <- approxfun(x = boundary$Month, y = boundary$Uw, method = "linear", rule = 2)
vt <- approxfun(x = boundary$Month, y = boundary$vt, method = "linear", rule = 2)


times <- seq(from = 0, to = 12, by = 1)
yini <- c(5,5)

out <- ode(times = times, y = yini, func = TwoLayer, parms = parameters, method = 'rk4')

