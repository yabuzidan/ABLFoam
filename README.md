# ABLFoam
Repository for neutral Atmospheric Boundary Layer modeling in OpenFoam 7 (www.openfoam.org) and OpenFoam 2206

Example cases from Bellegoni et al 2023

## Contents
1. `kEpsilon/`: modeling of neutral ABL with k-Epsilon turbulence model \*,
2. `kOmegaSST/`: modeling of neutral ABL with k-OmegaSST turbulence model \**,
3. `kOmegaSST_F1blend/`: modeling of neutral ABL with k-OmegaSST turbulence model and F1-blending \**.

\* _More information available in_:  
```
@article{PARENTE2011,
title = {Comprehensive Modelling Approach for the Neutral Atmospheric Boundary Layer: Consistent Inflow Conditions, Wall Function and Turbulence Model.},
journal = {Boundary-Layer Meteorology},
volume = {140},
pages = {411–428},
year = {2011},
doi = {https://doi.org/10.1007/s10546-011-9621-5},
author = {Alessandro Parente and Catherine Gorlé and Jeroen van Beeck and Carlo Benocci},
}
```
```
@article{LONGO2017160,
title = {Advanced turbulence models and boundary conditions for flows around different configurations of ground-mounted buildings},
journal = {Journal of Wind Engineering and Industrial Aerodynamics},
volume = {167},
pages = {160-182},
year = {2017},
doi = {https://doi.org/10.1016/j.jweia.2017.04.015},
author = {Riccardo Longo and Marco Ferrarotti and Clara García Sánchez and Marco Derudi and Alessandro Parente}
}
```
    
\** _More information available in_:  
```   
@article{BELLEGONI2023105583,
title = {An extended SST k−ω framework for the RANS simulation of the neutral Atmospheric Boundary Layer},
journal = {Environmental Modelling & Software},
volume = {160},
pages = {105583},
year = {2023},
doi = {https://doi.org/10.1016/j.envsoft.2022.105583},
author = {Marco Bellegoni and Léo Cotteleer and Sampath Kumar {Raghunathan Srikumar} and Gabriele Mosca and Alessandro Gambale and Leonardo Tognotti and Chiara Galletti and Alessandro Parente}
}
```
