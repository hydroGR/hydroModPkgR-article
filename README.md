# hydroModPkgR-article

Example codes to run 8 [hydrological modelling R packages](https://CRAN.R-project.org/view=Hydrology).
Each script enables application of the models on a simple hydrology example.

The input datasets are not provided, but each package includes at least one example dataset.

Only one parameter set is tested. Parameter estimation and uncertainty analysis procedures can be applied within the R environnement using the available features (embedded in some of the hydrological modelling packages or in other packages).

All scripts are based on the same pattern:

1. data formatting
2. parameter set definition
3. evaluation function definition
4. model simulation run
5. model run evaluation


Package      | Package references |
|---         |---         |
airGR        | [Coron et al., 2017](https://doi.org/10.1016/j.envsoft.2017.05.002); [Coron et al., 2020](https://doi.org/10.15454/EX11NA) |
dynatopmodel | [Metcalfe et al., 2015](https://doi.org/10.1016/j.envsoft.2015.06.010); [Metcalfe et al., 2018](https://cran.r-project.org/package=dynatopmodel)|
HBV.IANIGLA  | [Toum, 2019](https://CRAN.R-project.org/package=HBV.IANIGLA)
hydromad     | [Andrews et al., 2011](https://doi.org/10.1016/j.envsoft.2011.04.006); [Andrews and Guillaume, 2018](http://hydromad.catchment.org/)|
sacsmaR      | [Taner, 2019](https://github.com/tanerumit/sacsmaR)|
topmodel     | [Buytaert, 2018](https://CRAN.R-project.org/package=topmodel)|
TUWmodel     | [Parajka et al., 2007](https://doi.org/10.1002/hyp.6253); [Viglione and Parajka, 2020](https://CRAN.R-project.org/package=TUWmodel) |
WALRUS       | [Brauer et al., 2014a](https://doi.org/10.5194/gmd-7-2313-2014), [2014b](https://doi.org/10.5194/hess-18-4007-2014); [Brauer et al., 2017](https://github.com/ClaudiaBrauer/WALRUS) |

To cite these scripts: https://doi.org/10.15454/3PPKCL 

For more details about hydrological modelling R packages, see [Astagneau et al., 2020](https://hess.copernicus.org/preprints/hess-2020-498/).

