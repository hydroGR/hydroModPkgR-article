## =========================================================================
##
## Description : Run the Dynamic TOPMODEL rainfall-runoff model of the 
##               dynatopmodel paxkage on a midland catchment
##
## Authors     : Paul Astagneau   <paul.astagneau@inrae.fr>
##               Olivier Delaigue <olivier.delaigue@inrae.fr>
##               Guillaume Thirel <guillaume.thirel@inrae.fr>
##
## Date        : 2020-09-20
##
## Package     : dynatopmodel v1.2.1
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
# install.packages("dynatopmodel")

## Loading packages
library(dynatopmodel)
library(raster)
library(xts)

## Set time zone because the dynatopmodel package was made with another 
## time reference
Sys.setenv(TZ = "UTC") 

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

## Digital elevation model
DEM <- raster("DEM_midland_25m.tif")





## =========================================================================
## -------- Step 1.1: preparing data format required by dynatopmodel
## =========================================================================

## Dates 
basinObsTS$formatDate <- as.POSIXct(as.character(basinObsTS$Date), "%Y%m%d", tz = "UTC")

## xts zoo objects and values in m/hr
PFordynatopmodel <- xts(x         = ((basinObsTS$Ptot)*0.001/24), 
                        order.by  = basinObsTS$formatDate, 
                        frequency = 1)

EFordynatopmodel <- xts(x         = ((basinObsTS$PET)*0.001/24), 
                        order.by  = basinObsTS$formatDate, 
                        frequency = 1) 

QFordynatopmodel <- xts(x         = ((basinObsTS$Qmm)*0.001/24), 
                        order.by  = basinObsTS$formatDate, 
                        frequency = 1)  





## =========================================================================
## -------- Step 1.2: time period
## =========================================================================

start <- "1990-01-01"
end   <- "2009-12-31"

## Selecting data on the specified time periods
P    <- window(PFordynatopmodel, start = start, end = end) 

PET  <- window(EFordynatopmodel, start = start, end = end) 

Qobs <- window(QFordynatopmodel, start = start, end = end) 

## Check data gaps 
print(summary(P))
print(summary(PET))
print(summary(Qobs))





## =========================================================================
## -------- Step 1.3: spatial distribution preprocessing work
## =========================================================================

## If the DEM includes sinks you can use the following function
## layers <- build_layers(DEM, fill.sinks = TRUE, deg = 40) 

## With a "clean" DEM
## The DEM has to have the same x and y resolution
## This function calculates the topographic index and 
## the upslope contributing areas for each cell of the DEM
## Warning: this function takes time to process
layers <- dynatopmodel::build_layers(DEM)

## This function identifies the cell of the DEM that contains the channel
## It is more accurate with a digital river network
## Without a DRN the cells containing the channel are the cells 
## with a topographic index greater than atb.thresh
chans <- dynatopmodel::build_chans(drn = NULL, atb = layers$atb, atb.thresh = 0.55) 

## Display the channel cells created by chans to check whether atb.thresh has a reasonable value
plot(chans)


## Here, the HRUs are created with the upslope contributing areas 
## With addLayers it could be other characteristics using raster layers such as flow distance layers
## Combinations can be made with the argument cuts=c() which defines the number of HRUs
## a is the name of the upslope contributing areas created by build_layers
disc <- dynatopmodel::discretise(layers, cuts = c(a = 5), chans = chans, area.thresh = 0.5/100) 


## Routing function
## The number of breaks can be changed if needed
## Warning: this function takes time to process
routing <- dynatopmodel::build_routing_table(dem = DEM, chans = chans$chans, breaks = 5) 





## =========================================================================
## -------- Step 2: parameters
## =========================================================================

## Several methods exist to perform deterministic or probabilistic parameter 
## estimation 
## Here one set of parameters is tested

## Parameters can be set differently for each HRU but will be set the same here
# vof       Overland flow velocity (m/h)
# m         Form of exponential decline in conductivity (m)
# srz_max   Max root zone storage (m)
# srz0      Initial root zone storage (-)
# vchan     Channel routing velocity (m/h)
# ln_t0     Lateral saturated transmissivity (m2/h)
# sd_max    Max effective deficit of saturated zone (m)
# td        Unsaturated zone time delay (h/m)

param <- c(vof   = 28,    m     = 0.025, srz_max = 0.12, srz0 = 0.98,
           vchan = 1755,  ln_t0 = 5.1,   sd_max  = 0.46, td   = 99)

## Assign the parameters to the HRUs
disc$groups$vof     <- param[1]
disc$groups$m       <- param[2]
disc$groups$srz_max <- param[3]
disc$groups$srz0    <- param[4]
disc$groups$vchan   <- param[5]
disc$groups$ln_t0   <- param[6]
disc$groups$sd_max  <- param[7]
disc$groups$td      <- param[8]





## =========================================================================
## -------- Step 3: create KGE function for evaluation of the outputs 
##                  (many other evaluations are possible such as uncertainty analysis)
## =========================================================================

## See (Gupta et al., 2009) for more precisions
kgeEval <- function(simu, obse, warmup) {
  ## warmup is a boolean for considering a 1 year warm-up period
  ## simu is the output of run.dtm (xts zoo object)
  ## obse is an xts zoo object
  
  obse <- as.numeric(obse[, 1])      # because runoff is an xts zoo object
  simu <- as.numeric(simu$qsim[, 1]) # because qsim is an xts zoo object
  simu <- simu[1:length(obse)]      # because the length of qsim sometimes differ from the length of Qobs
  
  if (warmup == TRUE) {
    simu <- simu[-c(1:365)]  # Removes the warm-up period of 1 year
    obse <- obse[-c(1:365)]  # Removes the warm-up period of 1 year
  } 
  
  # Calculations for KGE
  simu[is.na(obse)] <- NA
  mObse <- mean(obse, na.rm = TRUE)
  mSimu <- mean(simu, na.rm = TRUE)
  ere   <- sum((obse-mObse)*(simu-mSimu), na.rm = TRUE) / (sqrt(sum((obse-mObse)^2, na.rm = TRUE))*sqrt(sum((simu-mSimu)^2, na.rm = TRUE)))
  alpha <- sum((simu-mSimu)^2, na.rm = TRUE) / sum((obse-mObse)^2, na.rm = TRUE)
  beta  <- sum(simu, na.rm = TRUE) / sum(obse, na.rm = TRUE)
  sqrt(((ere-1)^2) + ((alpha-1)^2) + ((beta-1)^2)) 
  
}





## =========================================================================
## -------- Step 4: run Dynamic TOPMODEL with the selected set of parameters and visualize the results
## =========================================================================

## Run the model
runResults <- dynatopmodel::run.dtm(groups        = disc$groups, 
                                    weights       = disc$weights, 
                                    rain          = P, #rain in m/hr
                                    routing       = routing, 
                                    qobs          = Qobs, 
                                    qt0           = as.numeric(Qobs[1,]),
                                    pe            = PET,
                                    dt            = 24, 
                                    graphics.show = FALSE, 
                                    max.q         = 1 ) 

## KGE value
message(paste("KGE = ", round(1 - kgeEval(obse   = Qobs, 
                                          simu   = runResults, 
                                          warmup = TRUE),
                              digits = 2)))

## Plot the results
indicesPlot <- seq(from = which(format(basinObsTS$Date, format = "%Y%m%d") == gsub("-", "", start)), 
                   to  = which(format(basinObsTS$Date, format = "%Y%m%d") == gsub("-", "", end)))[-c(1:365)] # Without the warm-up period

dynatopmodel::disp_output(qobs  = Qobs[-c(1:365)]*1000, # *1000 to plot the results in mm/hr
                          qsim  = runResults$qsim[1:length(Qobs[-c(1:365)]),1]*1000, 
                          rain  = P[-c(1:365)]*1000,
                          evap  = runResults$ae[1:length(Qobs[-c(1:365)]),1]*1000,
                          max.q = 2, cex.main = 1, col.axis = "slategrey", las.time = 1)


## Outputs
str(runResults)


