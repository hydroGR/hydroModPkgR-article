## =========================================================================
##
## Description : Run the SAC-SMA-SNOW17 rainfall-runoff model of the sacsmaR 
##               package on a mountainous catchment
##
## Authors     : Paul Astagneau   <paul.astagneau@inrae.fr>
##               Olivier Delaigue <olivier.delaigue@inrae.fr>
##               Guillaume Thirel <guillaume.thirel@inrae.fr>
##
## Date        : 2020-09-20
##
## Package     : sacsmaR v0.0.1
##
## Comments    :
##               Inputs  : P, PET, Temp, Qobs, average elevation
##               Outputs : Qsim 
##
## Reference   : DOI:
##
## =========================================================================







## =========================================================================
## -------- Step 0: preamble
## =========================================================================

## Install package
# install.packages("devtools")
# devtools::install_github("tanerumit/sacsmaR")

## Loading package
library(sacsmaR)

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

## Command line to find the julian day of the year required by sacsmaR 
basinObsTS$JD <- as.POSIXlt(as.character(basinObsTS$Date), tz = "UTC", format = "%Y%m%d")$yday + 1

basinObsTS$Date <- as.POSIXct(as.character(basinObsTS$Date), format = "%Y%m%d", tz = "UTC")



## Loading hypsometric data (quantiles of catchment elevation)
## basinHypsoQuantiles column name:
## $Elevation Elevation (m a.s.l.)
basinHypsoQuantiles <- read.table(file             = "elevation_quantiles.txt",
                                  header           = TRUE,
                                  sep              = "",
                                  dec              = ".",
                                  stringsAsFactors = FALSE)

## The average elevation is taken because 
## the snow accounting function is used in the lumped mode
averageElevation <- mean(basinHypsoQuantiles$Elevation) # The name of the hypsometric quantile column must be "Elevation"





## =========================================================================
## -------- Step 1.2: time period
## =========================================================================

start <- "19890101"
end   <- "20081230"

## ## To select the values corresponding to the chosen time period
indRun <- seq(from = which(format(basinObsTS$Date, format = "%Y%m%d") == start), 
              to   = which(format(basinObsTS$Date, format = "%Y%m%d") == end))

## Check data gaps 
print(summary(basinObsTS[indRun, ]))





## =========================================================================
## -------- Step 2: parameters
## =========================================================================

## Several methods exist to perform deterministic or probabilistic parameter 
## estimation 
## Here one set of parameters is tested

## Snow parameters of SNOW17
# SCF       Snow correction factor (-)                         
# PXTEMP    Temperature that separates rain from snow (°C)  
# MFMAX     Maximum melt factor (mm/°C/6 h)                 
# MFMIN     Minimum melt factor (mm/°C/6 h)                
# UADJ      Wind function factor (mm/mb) mb = millibar                  
# MBASE     Melt base temperature (°C)                      
# TIPM      Antecedent snow temperature index (-)              
# PLWHC     Liquid water holding capacity (%)               
# NMF       Maximum negative melt factor (mm/°C/6 h)        
# DAYGM     Average daily ground melt (mm/d)              


## Soil moisture accounting parameters
# UZTWM     Upper zone tension water capacity [mm]                              
# UZFWM     Upper zone free water capacity [mm]                                 
# LZTWM     Lower zone tension water capacity [mm]                              
# LZFPM     Lower zone primary free water capacity [mm]                         
# LZFSM     Lower zone supplementary free water capacity [mm]                   
# UZK       Upper zone free water lateral depletion rate [1/d]                  
# LZPK      Lower zone primary free water depletion rate [1/d]                 
# LZSK      Lower zone supplementary free water depletion rate [1/d]           
# ZPERC     Percolation demand scale parameter [-]                              
# REXP      Percolation demand shape parameter [-]                               
# PFREE     Percolating water split parameter (decimal fraction)                
# PCTIM     Impervious fraction of the watershed area (decimal fraction)        
# ADIMP     Additional impervious areas (decimal fraction)                      
# RIVA      Riparian vegetation area (decimal fraction)                              
# SIDE      The ratio of deep recharge to channel base flow [-]                      
# RSERV     Fraction of lower zone free water not transferrable (decimal fraction)  

param <- c(SCF   = 1.6,  PXTEMP = 1.1, MFMAX = 1.8,  MFMIN = 0.5,  UADJ  = 0.09,
           MBASE = 0.16, TIPM   = 0.9, PLWHC = 0.1,  NMF   = 0.1,  DAYGM = 0.02,
           UZTWM = 39,   UZFWM  = 12,  LZTWM = 11,   LZFPM = 194,  LZFSM = 902,  UZK  = 0.29, LZPK = 0.06, LZSK  = 0.01, 
           ZPERC = 90,   REXP   = 4,   PFREE = 0.68, PCTIM = 0.02, ADIMP = 0.07, RIVA = 0,    SIDE = 0.3,  RSERV = 0)  





## =========================================================================
## -------- Step 3: create KGE function for evaluation of the outputs 
##                  (many other evaluations are possible such as uncertainty analysis)
## =========================================================================

## See (Gupta et al., 2009) for more precisions
kgeEval <- function(obse, simu, warmup) {
  ## warmup is a boolean for considering a 1 year warm-up period
  ## simu is the output of sacSma (vector)
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
## -------- Step 4: run SAC-SMA with the selected set of parameters and visualize the results
## =========================================================================

## Run the model
precipLiquid <- sacsmaR::snow17(par   = param[1:10], # snow parameters
                                prcp  = basinObsTS$Ptot[indRun], 
                                tavg  = basinObsTS$Temp[indRun], 
                                elev  = averageElevation, 
                                jdate = basinObsTS$JD[indRun])

runResults <- sacsmaR::sacSma(par  = param[11:26], # hydro parameters
                              prcp = precipLiquid, 
                              pet  = basinObsTS$PET[indRun])$tot

## KGE value
message(paste("KGE = ", round(1-kgeEval(obse   = basinObsTS$Qmm[indRun], 
                                        simu   = runResults, 
                                        warmup = TRUE), 
                              digits = 2)))

## Plot results without warm-up period
## Observed hydrograph
plot(x = basinObsTS$Date[indRun[-c(1:365)]], 
     y = basinObsTS$Qmm[indRun[-c(1:365)]], 
     type = "l",
     xlab = "",
     ylab = "flow (mm/d)")

## Add simulated hydrograph
lines(x = basinObsTS$Date[indRun[-c(1:365)]], 
      y = runResults[-c(1:365)], 
      type = "l", 
      col = "red")

legend(x = "top", legend = c("obs", "sim"), col=c("black", "red"), lty = 1)


## Outputs
str(runResults)


