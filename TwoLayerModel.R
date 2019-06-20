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



TwoLayer <- function(t, y, bound, parms){
  dTe <-  Q / Ve * Tin -
    Q / Ve +
    ((vt * At) / Ve) * (dTh - dTe) +
    Jsn / (rho * cp * Ht) +
    (rho * (Tair + 273)^4 * (As + 0.031 * sqrt(eair)) * (1 - Rl)) / (rho * cp * Ht) -
    (eps * rho * (Te + 273)^4) /(rho * cp * Ht) -
    (c1 * Uw * (Te - Tair)) / (rho * cp * Ht) -
    (Uw * (es - eair)) / (rho * cp * Ht)
  dTh <-  (vt * At) / (Vh) * (Te - Th)
  
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
vt_mix <- 1

boundary <- cbind(boundary, c(rep(vt_mix,4),rep(vt,5),rep(vt_mix,3)))
colnames(boundary) <- c('Month','Jsw','Tair','Dew','U','vt')
