## =========================================================================
##
## Description : Run the IHACRES-CMD rainfall-runoff model of the hydromad 
##               package on a mountainous catchment
##
## Authors     : Paul Astagneau   <paul.astagneau@inrae.fr>
##               Olivier Delaigue <olivier.delaigue@inrae.fr>
##               Guillaume Thirel <guillaume.thirel@inrae.fr>
##
## Date        : 2020-09-20
##
## Package     : hydromad v0.9.26
##
## Comments    :
##               Inputs  : P, PET, Temp, Qobs
##               Outputs : Qsim 
##
## Reference   : DOI:
##
## =========================================================================







## =========================================================================
## -------- Step 0: preamble
## =========================================================================

## Install package
# install.packages(c("zoo", "latticeExtra", "polynom", "car", "Hmisc","reshape"))
# install.packages("hydromad", repos="http://hydromad.catchment.org")

## Loading packages
library(hydromad)
library(zoo)

## Loading hydroclimatic time series
## basinObsTS column names:
## $Date   YYYYMMDD
## $Ptot   Precipitations (liquid + solid) (mm/d)
## $PET    Potential evapotranspiration (mm/d)
## $Temp   Air temperature (°C)
## $Qmm    Runoff depth (mm/d) 
basinObsTS <- read.table(file = "mountainous_catchment_data_TS.txt",
                         header = TRUE,
                         sep = ";",
                         dec = ".")
basinObsTS$formatDate <- as.Date(as.character(basinObsTS$Date), "%Y%m%d", tz = "UTC") 





## =========================================================================
## -------- Step 1.1: preparing data format required by hydromad
## =========================================================================

PForHydromad <- zoo(x         = basinObsTS$Ptot, 
                    order.by  = basinObsTS$formatDate, 
                    frequency = 1)

EForHydromad <- zoo(x         = basinObsTS$PET, 
                    order.by  = basinObsTS$formatDate, 
                    frequency = 1) 

QForHydromad <- zoo(x         = basinObsTS$Qmm, 
                    order.by  = basinObsTS$formatDate, 
                    frequency = 1)

TForHydromad <- zoo(x         = basinObsTS$Temp, 
                    order.by  = basinObsTS$formatDate, 
                    frequency = 1)

## hydromad input is a zoo object with the column names P, E, Q, T 
dataForHydromad <- merge(P = PForHydromad, 
                         E = EForHydromad, 
                         Q = QForHydromad,
                         T = TForHydromad)





## =========================================================================
## -------- Step 1.2: time period
## =========================================================================

start <- "1989-01-01"
end   <- "2008-12-31"

## To select the values corresponding to the chosen time period
dataPeriod <- window(dataForHydromad, start = start, end = end) 

## Check data gaps 
print(summary(dataPeriod))





## =========================================================================
## -------- Step 2: parameters
## =========================================================================

## Several methods exist to perform deterministic or probabilistic parameter 
## estimation 
## Here one set of parameters is tested

## Tmax   temperature threshold for rain, all rain is liquid above this threshold (°C)
## Tmin   temperature threshold for rain, all rain is snow below this threshold   (°C)
## cr     correction factor for rainfall (-)
## cs     correction factor for snowfal  (-)
## kd     degree day factor for snowmelt (mm/°C/d)
## kf     degree day factor for freezing (mm/°C/d)
## rcap   retention parameter for liquid water capacity of snowpack (-)
## f      catchment moisture deficit stress threshold as a proportion of d (-)
## e      temperature to PET conversion factor (-)
## d      catchment moisture deficit threshold for producing flow (mm)
## tau_s  time constant for the exponential component, slow response (d)
## tau_q  time constant for the exponential component, fast response (d)
## v_s    fractional volume for the exponential components, slow response (-)
param <- c(Tmax    = 0.262, Tmin    = 0.124, 
           cr      = 1.01,  cs      = 0.832, kd  = 2.06, kf   = 1.08,  rcap = 0.217,
           f       = 0.961, e       = 0.578, d   = 200,
           tau_s   = 30,    tau_q   = 5,     v_s = 0.5)





## =========================================================================
## -------- Step 3: create KGE function for evaluation of the outputs 
##                  (many other evaluations are possible such as uncertainty analysis)
## =========================================================================

## See (Gupta et al., 2009) for more precisions
hydromad::hydromad.stats("KGE" = function(Q, X, ...) {
  
  obse <- Q    
  simu <- X 
  
  # Calculations of KGE
  simu[is.na(obse)] <- NA
  mObse <- mean(obse, na.rm = TRUE)
  mSimu <- mean(simu, na.rm = TRUE)
  ere   <- sum((obse-mObse)*(simu-mSimu), na.rm = TRUE) / (sqrt(sum((obse-mObse)^2, na.rm = TRUE))*sqrt(sum((simu-mSimu)^2, na.rm = TRUE)))
  alpha <- sum((simu-mSimu)^2, na.rm = TRUE) / sum((obse-mObse)^2, na.rm = TRUE)
  beta  <- sum(simu, na.rm = TRUE) / sum(obse, na.rm = TRUE)
  1 - sqrt(((ere-1)^2) + ((alpha-1)^2) + ((beta-1)^2)) 
})





## =========================================================================
## -------- Step 4: run IHACRES-CMD with the selected set of parameters and visualize the results
## =========================================================================

## Several methods exist to perform deterministic or probabilistic parameter 
## estimation 
## Here one set of parameters is tested

## Run the model with the selected parameter set
runResults <- hydromad::hydromad(DATA    = dataPeriod, 
                                 sma     = "snow", # to use the snow routine with the CMD soil moisture accounting routine
                                 routing = "expuh",
                                 warmup  = 365,
                                 Tmax    = param[1],  Tmin    = param[2], 
                                 cr      = param[3],  cs      = param[4],  kd  = param[5], kf   = param[6],  rcap = param[7],
                                 f       = param[8],  e       = param[9],  d   = param[10],
                                 tau_s   = param[11], tau_q   = param[12], v_s = param[13]) 

## KGE value
message(paste("KGE =", round(hydromad::hmadstat("KGE")(Q = dataPeriod$Q[-c(1:365)],              # Without the warm-up period
                                                       X = runResults$fitted.values[-c(1:365)]), # Without the warm-up period
                             digits = 2))) 


## Plot the results
c(xyplot(runResults, ylab = "mm/d", with.P = TRUE),
  residuals = xyplot(residuals(runResults, type = "h"), ylab = "mm/d"),
  layout    = c(1,3), 
  y.same    = FALSE)


## Outputs
str(runResults$U)

