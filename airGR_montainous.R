## =========================================================================
##
## Description : Run the GR4J rainfall-runoff model with the CemaNeige snow 
##               model of the airGR package on a montainous catchment
##
## Authors     : Paul Astagneau   <paul.astagneau@inrae.fr>
##               Olivier Delaigue <olivier.delaigue@inrae.fr>
##               Guillaume Thirel <guillaume.thirel@inrae.fr>
##
## Date        : 2020-09-20
##
## Package     : airGR v1.4.3.65
##
## Comments    : 
##               Inputs  : P, PET, Temp, Qobs, hypsometric curve
##               Outputs : Qsim 
##
## Reference   : DOI:
##
## =========================================================================







## =========================================================================
## -------- Step 0: preamble
## =========================================================================

## Install package
# install.packages("airGR")

## Loading package
library(airGR)

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

## Loading hypsometric data (quantiles of catchment elevation)
## basinHypsoQuantiles column name:
## $Elevation
basinHypsoQuantiles <- read.table(file             = "elevation_quantiles.txt",
                                  header           = TRUE,
                                  sep              = "",
                                  dec              = ".",
                                  stringsAsFactors = FALSE)





## =========================================================================
## -------- Step 1.1: creating the inputs object for GR4J-CemaNeige
## =========================================================================

inputsModel <- airGR::CreateInputsModel(FUN_MOD   = RunModel_CemaNeigeGR4J, 
                                        DatesR    = basinObsTS$Date,
                                        Precip    = basinObsTS$Ptot, 
                                        PotEvap   = basinObsTS$PET,
                                        TempMean  = basinObsTS$Temp, 
                                        HypsoData = basinHypsoQuantiles$Elevation,
                                        ZInputs   = median(basinHypsoQuantiles$Elevation))





## =========================================================================
## -------- Step 1.2: time period
## =========================================================================

## airGR automatically deals with the warm-up period (default = 1 year) if 
## not user-defined
start <- "19900101"
end   <- "20081231"

## To select the indices corresponding to the chosen time period
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

## X1	production store capacity (mm)
## X2	intercatchment exchange coefficient (mm/d)
## X3	routing store capacity (mm)
## X4	unit hydrograph time constant (d)
## CNX1	weighting coefficient for snow pack thermal state (-)
## CNX2	degree-day melt coefficient (mm/°C/d)
param <- c(X1 = 393, X2 = 0.6, X3 = 143, X4 = 1.2, CNX1 = 0.4, CNX2 = 5.2)





## =========================================================================
## -------- Step 3: KGE function for evaluation of the outputs 
##                  (many other evaluations are possible such as uncertainty analysis)
## =========================================================================

## Prepare running options
runOptionsEval <- airGR::CreateRunOptions(FUN_MOD          = RunModel_CemaNeigeGR4J,
                                          InputsModel      = inputsModel, 
                                          IndPeriod_Run    = indRun,
                                          IniStates        = NULL, 
                                          IniResLevels     = NULL, 
                                          IndPeriod_WarmUp = NULL)

## Prepare KGE function 
## See (Gupta et al., 2009) for more precisions
inputsCritEval <- airGR::CreateInputsCrit(FUN_CRIT    = ErrorCrit_KGE, 
                                          InputsModel = inputsModel, 
                                          RunOptions  = runOptionsEval, 
                                          Obs         = basinObsTS$Qmm[indRun])





## =========================================================================
## -------- Step 4: run GR4J-CemaNeige with the selected set of parameters and visualize the results
## =========================================================================

## Run the model
runResults <- airGR::RunModel_CemaNeigeGR4J(InputsModel   = inputsModel,
                                            RunOptions  = runOptionsEval,
                                            Param       = param)





## =========================================================================
## -------- Step 5: model run evaluation
## =========================================================================

## Evaluate simulations by calculating the KGE criterion (many other evaluations are possible such as uncertainty analysis)
outputsCrit <- airGR::ErrorCrit_KGE(InputsCrit   = inputsCritEval,
                                    OutputsModel = runResults)


## Plot results 
plot(runResults, Qobs = basinObsTS$Qmm[indRun], which = "all")


## Outputs
str(runResults)
