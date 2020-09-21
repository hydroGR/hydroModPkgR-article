## =========================================================================
##
## Description : Run the WALRUS rainfall-runoff model of the WALRUS package 
##               on a midland catchment
##
## Authors     : Paul Astagneau   <paul.astagneau@inrae.fr>
##               Olivier Delaigue <olivier.delaigue@inrae.fr>
##               Guillaume Thirel <guillaume.thirel@inrae.fr>
##
## Date        : 2020-09-20
##
## Package     : WALRUS v1.10
##
## Comments    :
##               Inputs  : P, PET, Qobs
##               Outputs : Qsim
##
## Reference   : DOI:
##
## =========================================================================







## =========================================================================
## -------- Step 0: preamble
## =========================================================================

## Install package
# library(devtools)
# install_github("ClaudiaBrauer/WALRUS")

## Loading package
library(WALRUS)

## Loading hydroclimatic time series
## basinObsTS column names:
## $Date   YYYYMMDD
## $Ptot   Precipitations (liquid + solid) (mm/d)
## $PET    Potential evapotranspiration (mm/d)
## $Temp   Air temperature (°C)
## $Qmm    Runoff depth (mm/d) 
basinObsTS <- read.table(file = "midland_catchment_data_TS.txt",
                         header = TRUE,
                         sep = ";",
                         dec = ".")





## =========================================================================
## -------- Step 1.1: preparing data format required by WALRUS
## =========================================================================

## The column names must be as below
## Data are in mm/timestep according to section 5.5 page 20 of the user manual
dataForWALRUS <- data.frame(date  = basinObsTS$Date,  P = basinObsTS$Ptot,
                            ETpot = basinObsTS$PET, Q = basinObsTS$Qmm)





## =========================================================================
## -------- Step 1.2: time period
## =========================================================================

## To select the values corresponding to the time dataPeriod
dataPeriod <- WALRUS::WALRUS_selectdates("dataForWALRUS", 19900100, 20091231)

## Define a 1 year warm-up dataPeriod
dataPeriod$warm <- 24*365 #in hours

## Check data gaps 
print(summary(dataPeriod))





## =========================================================================
## -------- Step 2: parameters
## =========================================================================

## Several methods exist to perform deterministic or probabilistic parameter 
## estimation 
## Here one set of parameters is tested

## cW wetness index parameter        (mm)
## cV vadose zone relaxation time    (h)
## cG groundwater reservoir constant (mm.h)
## cQ quickflow reservoir constant   (h)
## cS surface water parameter        (mm/h)
## CD average channel depth/bottom   (mm)
## as surface water area fraction (= fraction of catchment covered by ditches and channels) (-) 
## st soil type (-)
## CD, as and st are determined from catchment characteristics (see section 6.3 page 25 of the user manual)
param <- data.frame(cW = 315, cV = 0.1,  cG = 1e8,  cQ = 67, 
                    cS = 5,   cD = 2500, aS = 0.01, st = "loamy_sand", stringsAsFactors = FALSE)





## =========================================================================
## -------- Step 3: create KGE function for evaluation of the outputs 
##                  (many other evaluations are possible such as uncertainty analysis)
## =========================================================================

## See (Gupta et al., 2009) for more precisions
kgeEval <- function(obse, simu, warmup) {
  ## warmup is a boolean for considering a 1 year warm-up dataPeriod
  ## simu is the output of WALRUS_loop()$Q (vector)
  ## obse is a vector
  
  if (warmup == TRUE) {
    simu <- simu[-c(1:365)]          # Removes the warm-up dataPeriod of 1 year
    simu <- simu[1:(length(simu)-1)] # because of a not wanted additional value 
    obse <- obse[-c(1:365)]          # Removes the warm-up dataPeriod of 1 year
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
## -------- Step 4: run WALRUS with the selected set of parameters and visualize the results
## =========================================================================

WALRUS::WALRUS_preprocessing(f = dataPeriod, dt = 1) 
name <-  "one_parameter_set"

## Run the model 
runResults <- WALRUS::WALRUS_loop(pars = param)

## KGE value
message(paste("KGE =", round(1 - kgeEval(obse   = dataPeriod$Q,
                                         simu   = runResults$Q,
                                         warmup = TRUE), 
                             digits = 2)))

## Plot results without warm-up dataPeriod
## WARNING: There must be two folders, "figures" and "output" in the current directory
## This will also generate a plot in the R session
dir.create("figures")
dir.create("output")
WALRUS::WALRUS_postprocessing(o = runResults , pars = param , n = name, residuals = TRUE)


## Outputs
str(runResults)


