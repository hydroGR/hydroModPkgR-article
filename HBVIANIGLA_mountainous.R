## =========================================================================
##
## Description : Run the HBV rainfall-runoff model of the HBV.IANIGLA 
##               package on a mountainous catchment
##
## Authors     : Paul Astagneau   <paul.astagneau@inrae.fr>
##               Olivier Delaigue <olivier.delaigue@inrae.fr>
##               Guillaume Thirel <guillaume.thirel@inrae.fr>
##
## Date        : 2020-09-20
##
## Package     : HBV.IANIGLA v0.1.1
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
# install.packages("HBV.IANIGLA")

## Loading package
library(HBV.IANIGLA)

## Loading hydroclimatic time series
## basinObsTS column names:
## $Date   YYYYMMDD
## $Ptot   Precipitations (liquid + solid) (mm/d)
## $PET    Potential evapotranspiration (mm/d)
## $Temp   Air temperature (°C)
## $Qmm    Runoff depth (mm/d) 
basinObsTS <- read.table(file = "mountainous_catchment_data_TS.txt",
                         header = TRUE,
                         sep    = ";",
                         dec    = ".")
basinObsTS$Date <- as.POSIXct(as.character(basinObsTS$Date), format = "%Y%m%d", tz = "UTC")





## =========================================================================
## -------- Step 1.2: time period
## =========================================================================

start <- "19890101"
end   <- "20081231"

## To select the values corresponding to the chosen time period
indRun <- seq(from = which(format(basinObsTS$Date, format = "%Y%m%d") == start), 
              to   = which(format(basinObsTS$Date, format = "%Y%m%d") == end))

## Check data gaps 
print(summary(basinObsTS[indRun, ]))

## Select corresponding data
inputPeriod <- matrix(data = c(basinObsTS$Temp[indRun], basinObsTS$Ptot[indRun]), 
                      ncol = 2, dimnames = list(NULL, c("airT", "precip")))




## =========================================================================
## -------- Step 2: parameters
## =========================================================================

## Several methods exist to perform deterministic or probabilistic parameter 
## estimation 
## Here one set of parameters is tested

## FC    maximum soil moisture storage (mm)
## LP    parameter related to the limit for potential evaporation (-)
## Beta  the non linear parameter for runoff production (-)
## K0    storage coefficient for very fast response (1/d)
## K1    storage coefficient for fast response (1/d)
## K2    storage coefficient for slow response (1/d)
## UZL   threshold storage state (mm)
## PERC  constant percolation rate (mm/d)
## Bmax  base time of the unit hydrograph (d)
## SFCF  snowfall correction factor (-)
## Tr    solid and liquid precipitation threshold temperature (°C)
## Tt    melt temperature (°C)
## fm    snowmelt factor (mm/°C/d)
## UZL must always be greater than PERC and 1 > K0 > K1 > K2 
param <- c(FC  = 291, LP   = 0.8, Beta = 1.4, K0   = 0.5,  K1 = 0.13, K2 = 0.03, 
           UZL = 62,  PERC = 4.1, Bmax = 1.7, SFCF = 0.97, Tr = 1.3,  Tt = -1.3, fm = 1)





## =========================================================================
## -------- Step 3: create KGE function for evaluation of the outputs 
##                  (many other evaluations are possible such as uncertainty analysis)
## =========================================================================

## See (Gupta et al., 2009) for more precisions
kgeEval <- function(obse, simu, warmup) {
  ## warmup is a boolean for considering a 1 year warm-up period
  ## simu is the output of UH (vector)
  ## obse is a vector
  
  if (warmup == TRUE) {
    simu <- simu[-c(1:365)]  # Removes the warm-up period of 1 year
    obse <- obse[-c(1:365)] # Removes the warm-up period of 1 year
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
## -------- Step 4: run HBV with the selected set of parameters and visualize the results
## =========================================================================

## Run the model 
resSnow <- HBV.IANIGLA::SnowGlacier_HBV(model     = 1,  # To select the HBV temperature index model
                                        inputData = inputPeriod, 
                                        initCond  = c(10, 2, 0) , 
                                        param     = param[10:13])

resProd <- HBV.IANIGLA::Soil_HBV(model     = 1,         # To select the classical soil moisture HBV routine
                                 inputData = matrix(data     = c(resSnow[, "Total"], basinObsTS$PET[indRun]), 
                                                    ncol     = 2, 
                                                    dimnames = list(NULL, c("Total", "PET"))), 
                                 initCond  = c(50, 1), 
                                 param     = param[1:3])

resRout <- HBV.IANIGLA::Routing_HBV(model     = 3,      # To select the HBV routing routine with two reservoirs and three flow outlets
                                    lake      = FALSE, 
                                    inputData = matrix(data = resProd[, "Rech"], ncol = 1), 
                                    initCond  = c(10, 15), 
                                    param     = param[4:8])

runResults <- HBV.IANIGLA::UH(model = 1,                # To select the triangular delay function
                              Qg    = matrix(data = resRout[, "Qg"], ncol = 1), 
                              param = param[9])

## KGE value
message(paste("KGE = ", round(1-kgeEval(simu   = runResults,
                                        obse   = basinObsTS$Qmm[indRun], 
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
      y    = runResults[-c(1:365)], 
      type = "l", 
      col  = "red")

legend(x = "top", legend = c("obs", "sim"), col = c("black", "red"), lty = 1)


## Outputs
str(resSnow)
str(resProd)
str(resRout)



