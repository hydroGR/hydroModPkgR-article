##  ========================================================================= 
##
## Description : Run the (lumped) HBV rainfall-runoff model of the TUWmodel 
##               package on a montainous catchment
##
## Authors     : Paul Astagneau   <paul.astagneau@inrae.fr>
##               Olivier Delaigue <olivier.delaigue@inrae.fr>
##               Guillaume Thirel <guillaume.thirel@inrae.fr>
##
## Date        : 2020-09-20
##
## Package     : TUWmodel v1.1-1
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
# install.packages("TUWmodel")

## Loading package
library(TUWmodel)

## Loading hydroclimatic time series
## basinObsTS column names:
## $Date   YYYYMMDD
## $Ptot   Precipitations (liquid + solid) (mm/d)
## $PET    Potential evapotranspiration (mm/d)
## $Temp   Air temperature (°C)
## $Qmm    Runoff depth (mm/d) 
basinObsTS <- read.table(file   = "mountainous_catchment_data_TS.txt", 
                         header = TRUE, 
                         sep    = ";", 
                         dec    = ".")
basinObsTS$Date <- as.POSIXct(as.character(basinObsTS$Date), format = "%Y%m%d", tz = "UTC")





##  ========================================================================= 
## -------- Step 1.2: time period
##  ========================================================================= 

start <- "19890101"
end   <- "20081231"

## To select the values corresponding to the time period
indRun <- seq(from = which(format(basinObsTS$Date, format = "%Y%m%d") == start), 
              to   = which(format(basinObsTS$Date, format = "%Y%m%d") == end))

## Check data gaps 
print(summary(basinObsTS[indRun, ]))





##  ========================================================================= 
## -------- Step 2: parameters
##  ========================================================================= 

## Several methods exist to perform deterministic or probabilistic parameter 
## estimation 
## Here one set of parameters is tested

## SCF    snow correction factor (-) 
## DDF    degree day factor (mm/degC/day) 
## Tr     threshold temperature above which precipitation is rain (degC) 
## Ts     threshold temperature below which precipitation is snow (degC) 
## Tm     threshold temperature above which melt starts (degC) 
## LPrat  parameter related to the limit for potential evaporation (-) 
## FC     field capacity, i.e., max soil moisture storage (mm) 
## BETA   the non linear parameter for runoff production  (-) 
## k0     storage coefficient for very fast response      (days) 
## k1     storage coefficient for fast response (days) 
## k2     storage coefficient for slow response (days) 
## lsuz   threshold storage state, i.e., the very fast response start if exceeded (mm) 
## cperc  constant percolation rate (mm/day) 
## bmax   maximum base at low flows (days) 
## croute free scaling parameter    (days^2/mm)
param <- c(SCF  = 1.2,  DDF = 0.9, Tr = 2.5, Ts = -2.3, Tm   = -1.6, LPrat = 0.2, FC = 190,
           BETA = 3.5,  k0  = 0.6, k1 = 12,  k2 = 38,   lsuz = 93,   cperc = 3.8, bmax = 2.2, croute = 0.6)





## =========================================================================
## -------- Step 3: create KGE function for evaluation of the outputs 
##                  (many other evaluations are possible such as uncertainty analysis)
## =========================================================================

## See (Gupta et al., 2009) for more precisions
kgeEval <- function(obse, simu, warmup) {
  ## warmup is a boolean for considering a 1 year warm-up period
  ## simu is the output of TUWmodel()$q[1,] (vector)
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





##  ========================================================================= 
## -------- Step 4: run TUWmodel with the selected set of parameters and visualize the results
##  ========================================================================= 

## Run the model 
runResults <- TUWmodel::TUWmodel(prec  = basinObsTS$Ptot[indRun], 
                                 airt  = basinObsTS$Temp[indRun], 
                                 ep    = basinObsTS$PET[indRun],
                                 area  = 1, # to select the lumped model
                                 param = param)





## =========================================================================
## -------- Step 5: model run evaluation
## =========================================================================

## KGE value
message(paste("KGE = ", round(1 - kgeEval(obse   = basinObsTS$Qmm[indRun],
                                          simu   = runResults$q[1,],
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
      y    = runResults$q[1, -c(1:365)], 
      type = "l", 
      col  = "red")

legend(x = "top", legend = c("obs", "sim"), col = c("black", "red"), lty = 1)


## Outputs
str(runResults)


