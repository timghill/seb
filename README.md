# SEB-clean
This repository contains surface energy balance (SEB) and subsurface model (SSM) code written in Matlab, developed for Kaskawulsh and Naluday (Lowell) Glaciers, St. Elias Mountains, Yukon.

The code is in the `code/` directory, while the `examples/` directory includes the necessary input files and script to model melt on Kaskawulsh Glacier in the 2018 melt season.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3923034.svg)](https://doi.org/10.5281/zenodo.3923034)


## Installation
Simply clone the repository and add the path to your Matlab path,
```
addpath(genpath('/path/to/repo/'))
```

## Running the model
The model is run using the `run_seb.m` script in the `code/` directory. This script takes as input a single structure `sim`, with fields given in Table 1. Table 2 describes the meteorological forcing the model requires to run, and table 3 lists all default parameters and options.

The model saves outputs with a frequency given by `sim.params.deltat`. Each output has fields given in Table 4. The model outputs melt in units of meters of ablation, **not** meters of water equivalent. To convert outputs to meters of water equivalent, multiple by `sim.params.rhoice/1000`.

**Table 1**: Attributes of `sim` structure to run the SEB model.

Attribute | Description
--------- | ------------------
`sim.paths.output` | Directory to save outputs into. Directory must exist.
`sim.paths.dem` | Path to `.mat` file containing digital elevation model. Must have attributes `dem.Z`, surface elevation in m; `dem.l`, DEM grid spacing in m; `dem.xx` and `dem.yy`, arrays of x and y coordinates in m.
`sim.paths.albedo` | Optional: Path to `.mat` file containing surface albedo map (must be same shape as DEM). If not given, model uses constant albedo given by `sim.params.albedo` (see default values table)
`sim.forcing` | Matlab structure with attributes given in Table 2.
`sim.params` | Structure containing any non-default parameters to run model with (see default values table). You will typically specify at least `sim.params.lat`, the latitude of the glacier; `sim.params.tout`, the frequency of model outputs; and `sim.params.deltat`, the frequency of the input meteorological forcing.

**Table 2**: Meteorological attributes required to run the model, including units.

Attribute | Units | Description
----------|-------|--------------------
`forcing.T` | <sup>o</sup>C (celsius) | Near-surface air temperature
`forcing.RH` | % [0-100] | Near-surface relative humidity
`forcing.LWin` | W m<sup>-2</sup> | Incoming Longwave radiation at the elevation of the temperature and RH measurements
`forcing.SWin` | W m<sup>-2</sup> | Incoming shortwave radiation incident on a horizontal surface
`forcing.P` | Pa | Air pressure at elevation of temperature measurements and with the specified temperature
`forcing.u` | m s<sup>-1</sup> | Average wind speed
`forcing.tt` | --| array of `datetime` objects representing the times of observations
`forcing.t` | hours | array of the decimal local solar time. E.g. 2:30 pm = 14.5.

**Table 3**: Default values and parameters for `sim.params`

Attribute | Units | Default value | Description
----------|-------|---------------|-----------------
`sigma` | W m<sup>-2</sup>K<sup>-4</sup> | 5.670374e-8 | Boltzmann constant
`eps_s` | -- | 1.0 | Ice surface emissivity
`cp` | J kg<sup>-1</sup> K<sup>-1</sup> | 1.005e3 | Specific heat capacity of air
`cp_ice` | J kg<sup>-1</sup> K<sup>-1</sup> | 2.108e3 | Specific heat capacity of ice
`cond` | W m<sup>-1</sup> K<sup>-1</sup> | 2.09 | Thermal conductivity of ice
`diff_ice` | m<sup>2</sup> s<sup>-1</sup> | `cond/(rhoice*cp_ice)` | Thermal diffusivity of ice
`rhoice` | kg m<sup>-3</sup> | 850 | Density of surface layer of ice
`Lv` | J kg<sup>-1</sup> | 2.264705e6 | Latent heat of vaporization of water
`Lf` | J kg<sup>-1</sup> | Latent heat of fusion of water
`kvk`| -- | 0.4 | von Karman's constant
`Rd` | J kg<sup>-1</sup> K<sup>-1</sup> | 287.058 | Dry air ideal gas constant
`Rv` | J kg<sup>-1</sup> K<sup>-1</sup> | 461.495 | Water vapour ideal gas constant
`g` | m s<sup>-2</sup> | 9.81 | Gravitational acceleration
`alpha` | -- | 0.35 | Ice surface albedo, if not using spatially distributed albedo as an input
`delta` | m | 5 | Constant used in surface shading algorithm to remove grid artifacts
`T_lapse_rate` | <sup>o</sup>C m<sup>-1</sup> | 3.9770e-3 | Temperature lapse rate (Positive means temperatures decrease with elevation)
`z0` | m | 0.003 | Momentum roughness length
`z0H` | m | 0.00003 | Heat roughness length
`z0E` | m | 0.00003 | Moisture roughness length
`z` | m | 1.5 | Height above ice surface of temperature and RH measurements
`T_elev` | m (asl.) | 760.153 | Elevation of temperature and RH measurements
`H` | m | 12 | Depth below surface of deep ice temperature ` Tice` used as boundary condition in subsurface model
`Tice` | <sup>o</sup>C | -3 | Ice temperature at depth `H`
`N_ssf` | -- | 12 | Number of vertical layers in subsurface model
`dt_ssm` | s | 900 | Time step for subsurface model
`cast_shadows` | `bool` | `true` | Switch to compute shading of glacier surface
`run_subsurface_model` | `bool` | `true` | Switch to run subsurface model

**Table 4**: Attributes of model outputs.

Attribute | Units |  Description
----------|-------|-----------------
`melt` | m | Cumulative surface ablation
`melt_SW` | m |Cumulative surface ablation due to shortwave radiation
`melt_LW` | m |Cumulative surface ablation due to longwave radiation
`melt_QE` | m |Cumulative surface ablation due to sensible heat flux
`melt_SH` | m |Cumulative surface ablation due to latent heat flux
`time` | -- | `datetime` object
`Qnet` | W m<sup>-2</sup> | Net energy balance
`LWnet` | W m<sup>-2</sup> | Net longwave radiation
`SWnet` | W m<sup>-2</sup> | Net shortwave radiation
`QE` | W m<sup>-2</sup> | Net sensible heat flux
`QH` | W m<sup>-2</sup> | Net latent heat flux
`SWin` | W m<sup>-2</sup> | Direct incident shortwave radiation
`LWin` | W m<sup>-2</sup> | Direct incident longwave radiation
`LWout` | W m<sup>-2</sup> | Outgoing longwave radiation
`Ts` | <sup>o</sup>C | Mean surface temperature
`albedo` | -- | Surface albedo
`QT` | W m<sup>-2</sup> | Warming heat flux
`QM` | W m<sup>-2</sup> | Melting heat flux

## Examples
The `examples/kaskawulsh/` directory contains an example case and is a good place to start. This case models melt on Kaskawulsh Glacier in 2018, not including surface shading since that makes the model run much slower.
