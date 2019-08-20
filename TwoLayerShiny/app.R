library(shiny)
library(deSolve)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)

ui <- fluidPage(
  titlePanel("Two-Layer Lake Model (Chapra 1997)"),
  'Model to simulate thermal regimes in lakes using a two layer assumption. Equations and derivations from Steven C. Chapra (2008): Surface Water-Quality Modeling. Waveland Press, Inc.',
  strong("Move the sliders to add/subtract values from the daily driver data:"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("solarInput", "Solar radiation", min = -100, max = 500, value = 0, post = " cal/(cm2 d)"),
      sliderInput("tempInput", "Air temperature", min = -5, max = 10, value = 0, post = " deg C"),
      sliderInput("windInput", "Wind speed", min = -100, max = 100, value = 0, post = " km/h"),
      sliderInput("thermoInput", "Thermocline transfer",min =-5, max = 50, value = 0, post = " m/d"),
      sliderInput("inflowInput", "Inflow temperature",min =-10, max = 25, value = 0, post = " deg C"),
      sliderInput("volumeInput", "Bottom Layer volume",min =5000, max = 175000, value = 75000, post = " m3")
    ),
    mainPanel(
      plotOutput("lakeplot")
    )
  )
)

server <- function(input, output) {
  output$lakeplot <- renderPlot({
    
    TwoLayer <- function(t, y, parms){
      eair <- (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t)))) # air vapor pressure
      esat <- 4.596 * exp((17.27 * Tair(t)) / (237.3 + Tair(t))) # saturation vapor pressure
      RH <- eair/esat *100 # relative humidity
      es <- 4.596 * exp((17.27 * y[1])/ (273.3+y[1]))
      # epilimnion water temperature change per time unit
      dTe <-  Q / Ve * Tin -              # inflow heat
        Q / Ve * y[1] +                   # outflow heat
        ((vt(t) * At) / Ve) * (y[2] - y[1]) + # mixing between epilimnion and hypolimnion
        + As/(Ve * rho * cp) * (
          Jsw(t)  + # shortwave radiation
            (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
            (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
            (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
            (Uw(t) * ((es) - (eair))) )# evaporation
      # hypolimnion water temperature change per time unit
      dTh <-  ((vt(t) * At) / Vh) * (y[1] - y[2]) 
      return(list(c(dTe, dTh)))
    }
    
    # input data: month, shortwave radiation, air temperature, dew point temperature, wind speed
    bound <- matrix(c(seq(1,12,1),
                      169, 274, 414, 552, 651, 684, 642, 537, 397, 259, 160, 127,
                      8.3, 9., 13.5,13.9,21.8,24.7,29.4,26.6,24.9,15.,9.7,6.6,
                      2.8,3.3,4.9,4.,5.3,7.8,11.8,11.5,7.7,6.8,6.5,2.4,
                      11.6,11.7,16.4,15.6,16.6,16.7,12.7,11.7,14.,12.9,14.8,11.6), nrow = 12, byrow = FALSE)
    Ve <- 175000 *1e6 # epilimnion volume
    Vh <- input$volumeInput * 1e6#75000*1e6 # hypolimnion volume
    At <- 11000 *1e4 # thermocline area
    Ht <- 3 * 100 # thermocline thickness
    As <- 25000 * 1e4 # surface area
    # H <- As/Ve#Ve/As *100
    Tin <- 10 # inflow water temperature
    Q <- 7500 *1e6 # inflow discharge
    Rl <- 0.03 # reflection coefficient (generally small, 0.03)
    Acoeff <- 0.6 # coefficient between 0.5 - 0.7
    sigma <- 11.7 * 10^(-8) # cal / (cm2 d K4) or: 4.9 * 10^(-3) # Stefan-Boltzmann constant in [J (m2 d K4)-1]
    eps <- 0.97 # emissivity of water
    rho <- 0.9982 # density (g per cm3)
    cp <- 0.99 # specific heat (cal per gram per deg C)
    c1 <- 0.47 # Bowen's coefficient
    
    Et <- 7.07 * 10^(-4)  * ((Ve+Vh)/As/100)^(1.1505) # vertifcal diffusion coefficient (cm2 per d)
    vto <- Et/(Ht/100) * (86400/10000) #*100 # heat change coefficient across thermocline during stratified season
    vt_mix <- 100 # high heat change coefficient during overturn and mixing events
    bound <- as.data.frame(cbind(bound, c(rep(vt_mix,3),rep(vto,6),rep(vt_mix,3))))
    
    colnames(bound) <- c('Month','Jsw','Tair','Dew','Uw','vt')
    bound$Uw <- 19.0 + 0.95 * (bound$Uw * 1000/3600)^2 # function to calculate wind shear stress (and transforming wind speed from km/h to m/s)
    
    bound.m <- bound
    bound.m $Jsw <- bound.m $Jsw + input$solarInput
    bound.m $Tair <- bound.m $Tair + input$tempInput
    bound.m $Uw <- bound.m $Uw + (19.0 + 0.95 * (input$windInput * 1000/3600)^2)
    bound.m $vt<- bound.m $vt+ input$thermoInput
    Tin<- Tin +input$inflowInput
    
    parameters <- c(Ve, Vh, At, Ht, As, Tin, Q, Rl, Acoeff, sigma, eps, rho, cp, c1)
  
    boundary <- rbind(bound.m [1,],bound.m )
    boundary$Month <- cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30,31))
  
  # approximating all boundary conditions 
    Jsw <- approxfun(x = boundary$Month, y = boundary$Jsw, method = "linear", rule = 2)
    Tair <- approxfun(x = boundary$Month, y = boundary$Tair, method = "linear", rule = 2)
    Dew <- approxfun(x = boundary$Month, y = boundary$Dew, method = "linear", rule = 2)
    Uw <- approxfun(x = boundary$Month, y = boundary$Uw, method = "linear", rule = 2)
    vt <- approxfun(x = boundary$Month, y = boundary$vt, method = "linear", rule = 2)
  
  #times <- seq(from = 1, to = nrow(boundary), by = 1/30)
    times <- seq(from = 1, to = max(boundary$Month), by = 1)
    yini <- c(5,5) # initial water temperatures
  
  # runge-kutta 4th order
    out <- ode(times = times, y = yini, func = TwoLayer, parms = parameters, method = 'rk4')
  
    result <- data.frame('Time' = out[,1],
                       'WT_epi' = out[,2], 'WT_hyp' = out[,3])
    ggplot(result) +
      geom_line(aes(x=Time, y=WT_epi, col='Surface Mixed Layer')) +
      geom_line(aes(x=(Time), y=WT_hyp, col='Bottom Layer')) +
      labs(x = 'Simulated Time in days', y = 'Water temperature in deg C')  +
      theme_bw()+
      scale_color_manual(values=c("blue", "red")) +
      ylim(0,30)+
      theme(legend.position="bottom")
  })
}

shinyApp(ui = ui, server = server)