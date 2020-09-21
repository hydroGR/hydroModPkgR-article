## =========================================================================
##
## Description : Run the TOPMODEL rainfall-runoff model of the topmodel 
##               package on a midland catchment
##
## Authors     : Paul Astagneau   <paul.astagneau@inrae.fr>
##               Olivier Delaigue <olivier.delaigue@inrae.fr>
##               Guillaume Thirel <guillaume.thirel@inrae.fr>
##
## Date        : 2020-09-20
##
## Package     : topmodel v0.7.3
##
## Comments    :
##               Inputs  : P, PET, Qobs, DEM
##               Outputs : Qsim 
##
## Reference   : DOI:
##
## =========================================================================







## =========================================================================
## -------- Step 0: preamble
## =========================================================================

## Install package
# install.packages("topmodel")

## Loading packages
library(topmodel)
library(raster)

## Loading hydroclimatic time series
## basinObsTS column names:
## $Date   YYYYMMDD
## $Ptot   Precipitations (liquid + solid) (mm/d)
## $PET    Potential evapotranspiration (mm/d)
## $Temp   Air temperature (°C)
## $Qmm    Runoff depth (mm/d) 
basinObsTS <- read.table(file   = "midland_catchment_data_TS.txt", 
                         header = TRUE, 
                         sep    = ";", 
                         dec    = ".")
basinObsTS$Date <- as.POSIXct(as.character(basinObsTS$Date), format = "%Y%m%d", tz = "UTC")

## Digital elevation model
DEM <- raster("DEM_midland_25m.tif")





## =========================================================================
## -------- Step 1.2: time period
## =========================================================================

start <- "19900101"
end   <- "20091231"

## For selecting data on period 1
indRun <- seq(from = which(format(basinObsTS$Date, format = "%Y%m%d") == start), 
              to   = which(format(basinObsTS$Date, format = "%Y%m%d") == end))

## Check data gaps 
print(summary(basinObsTS[indRun, ]))





## =========================================================================
## -------- Step 1.3: spatial distribution preprocessing work:
##                  topographic index and delay function
## =========================================================================

## Display DEM to check resolution
DEM 
image(DEM)

## DEM must be a matrix
DEM <- as.matrix(DEM)

# Remove the values outside the catchment:
DEM[DEM == -9999] <- NA

## To remove sinks from the DEM
DEM <- topmodel::sinkfill(DEM, res = 25, degree = 0.1)

## Calculation of the topographic index for each grid of the DEM
topidxDEM <- topmodel::topidx(as.matrix(DEM), resolution = 25)

## Classes of topographic index to prepare the topidx argument of the 
## topmodel() function 
topidxClasses <- topmodel::make.classes(topidxDEM$atb, n = 15) 

## Delay function
n <- 5 # five classes
delayFunction      <- topmodel::flowlength(DEM)*25 
delayFunction      <- topmodel::make.classes(delayFunction, n)
delayFunction      <- delayFunction[n:1,]
delayFunction[, 2] <- c(0, cumsum(delayFunction[1:(n-1), 2]))





## =========================================================================
## -------- Step 2: parameters
## =========================================================================

## Several methods exist to perform deterministic or probabilistic parameter 
## estimation 
## Here one set of parameters is tested

## qs0   	Initial subsurface flow per unit area (m)
## lnTe   log of the areal average of T0 (m2/h)
## m      Model parameter controlling the rate of decline of transmissivity in the soil profile (-)
## Sr0    Initial root zone storage deficit (m)
## SrMax  Maximum root zone storage deficit (m)
## td     Unsaturated zone time delay per unit storage deficit (h/m)
## vch    channel flow outside the catchment (m/h) 
## vr     channel flow inside catchment  (m/h)
## k0     Surface hydraulic conductivity (m/h)
## CD     capillary drive (-)
param <- c(qs0 = 0.00014, lnTe = 2.3,  m  = 0.024, Sr0 = 0.063, SrMax = 1.6, 
           td  = 0.76,    vch  = 1420, vr = 100,   k0  = 9,     CD    = 0.8)





## =========================================================================
## -------- Step 3: create KGE function for evaluation of the outputs 
##                  (many other evaluations are possible such as uncertainty analysis)
## =========================================================================

## See (Gupta et al., 2009) for more precisions
kgeEval <- function(obse, simu, warmup) {
  ## warmup is a boolean for considering a 1 year warm-up period
  ## simu is the output of topmodel (vector)
  ## obse is a vector
  
  if (warmup == TRUE) {
    simu <- simu[-c(1:365)]  # Removes the warm-up period of 1 year
    obse <- obse[-c(1:365)]  # Removes the warm-up period of 1 year
  }
  
  # Calculations for KGE
  simu[is.na(obse)] <- NA
  mObse <- mean(obse, na.rm = TRUE)
  mSimu <- mean(simu, na.rm = TRUE)
  ere   <- sum((obse-mObse)*(simu-mSimu), na.rm = TRUE) / (sqrt(sum((obse-mObse)^2, na.rm = TRUE))*sqrt(sum((simu-mSimu)^2, na.rm = TRUE)))
  alpha <- sum((simu-mSimu)^2, na.rm = TRUE)/sum((obse-mObse)^2, na.rm = TRUE)
  beta  <- sum(simu, na.rm = TRUE)/sum(obse, na.rm = TRUE)
  sqrt(((ere-1)^2) + ((alpha-1)^2) + ((beta-1)^2)) 
}





## =========================================================================
## -------- Step 4: run TOPMODEL with the selected set of parameters and visualize the results
## =========================================================================

## Run the model
param[11] <- 24  # time step expressed in hours
runResults <- topmodel::topmodel(parameters = param, 
                                 topidx     = topidxClasses, 
                                 delay      = delayFunction, 
                                 rain       = basinObsTS$Ptot[indRun]*0.001,   #in meters per time step
                                 ETp        = basinObsTS$PET[indRun]*0.001,    #in meters per time step
                                 verbose    = TRUE)  

## KGE value
message(paste("KGE = ", round(1-kgeEval(obse   = basinObsTS$Qmm[indRun],
                                        simu   = as.vector(runResults$Q)*1000, # in mm/d
                                        warmup = TRUE), 
                              digits = 2)))

## Plot results without warm-up period
## Observed hydrograph
plot(x    = basinObsTS$Date[indRun[-c(1:365)]], 
     y    = basinObsTS$Qmm[indRun[-c(1:365)]], 
     type = "l",
     xlab = "",
     ylab = "flow (mm/d)")

## Add simulated hydrograph
lines(x    = basinObsTS$Date[indRun[-c(1:365)]], 
      y    = as.vector(runResults$Q)[-c(1:365)]*1000, # in mm/day
      type = "l", 
      col  = "red")

legend(x = "top", legend = c("obs", "sim"), col=c("black", "red"), lty = 1)


## Outputs
str(runResults)






