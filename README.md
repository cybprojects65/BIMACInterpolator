# BIMAC - Bayesian Interpolation Model with Advection-diffusion Constraint

An advection equation solver for oceanic parameter spatial interpolation based on Jags and Gibbs sampling. Internally uses the stationary advection-diffusion equation as a constraint.

*How to cite:* G. Coro, (2023). An Open Science oriented Bayesian interpolation model for marine parameter observations. (to be published in) Environmental Modelling & Software.

![](https://github.com/cybprojects65/JagsOceanicSpatialInterpolator/blob/main/global_scale_example.png) *Figure 1. Example of global scale temperature layer estimated from observation data from the Argo network.*

The software is an all-in-one script containing all functions (BIMAC_functions.R). call_bimac_examples.R contains all examples. BIMAC comes with four internal models corresponding to different input types available to the user.

## Marine data interpolation with current velocity and depth constraints

*Required libraries:*

```         
  library (raster)
  library (R2jags)
  library (coda)
  library (plyr)
  library (dplyr)
  library (digest)
```

JAGS 4.3.1: <https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/>

*Input data:* Current-velocity u-v components' raster files in ASC format:

```         
currents_u_file = "./input_spatial_layers/oceanic_currents_u_wm.asc"

currents_v_file = "./input_spatial_layers/oceanic_currents_v_wm.asc" 
```

A file with the parameter observations to interpolate

```         
punctual_data_file = "./observations/temperature_argo.csv"
```

*Input parameters:* The bathymetry of the analysis (in meters)

```         
analysis_depth=-1
```

Number of adjacent cells to use for smoothing

```         
moving_average_points=1
```

An option to speed up the processing by using 100 iterations in the MCMC model.

```         
fast_solving=T
```

A parameter regulating the strictness to comply with the advection-diffusion equation: the smaller the stricter

```         
sd_advection_equation=0.1
```

Example:

```         
output<-bimac_currents_depthboundaries(punctual_data_file,
                                        currents_u_file,
                                        currents_v_file,
                                        analysis_depth=-1,
                                        moving_average_points=1, 
                                        fast_solving=T, 
                                        sd_advection_equation=0.1 )
file.rename("./output","./output_yes_currents_yes_depth")
```

The output is constituted of three files: An inverse-weighted distance interpolation serving as prior information for the MCMC model

```         
prior_data_output = "./output/BIMAC_IDW_prior.asc"
```

The final interpolation result (from the posterior distribution samples of the MCMC model)

```         
posterior_data_output = "./output/BIMAC_interpolation.asc"
```

The standard deviation associated to each interpolation cell

```         
posterior_data_output_sd = "./output/BIMAC_interpolation_sd.asc"
```

## Marine data interpolation with current velocity only

This function does not use bathymetry as a filter of land areas and interpolate over the entire bounding box enclosing the current-velocity files. The output is constituted of the same three files types described above.

Example:

```         
output<-bimac_nodepthboundaries(punctual_data_file,
                            currents_u_file,
                            currents_v_file, 
                            moving_average_points=1, 
                            fast_solving=T, 
                            sd_advection_equation=0.1 )
file.rename("./output","./output_yes_currents_no_depth")
```

## Marine data interpolation with depth constraint

This function only uses bathymetry as a constraint of land areas and does not require current-velocity files. The output is constituted of the same three files types described above.

Example:

```         
output<-bimac_noadvection(punctual_data_file,
                    analysis_depth=-1,
                    moving_average_points=1, 
                    fast_solving=T,
                    min_x_boundingbox=-0.25,
                    max_x_boundingbox=12,
                    min_y_boundingbox=34.75,
                    max_y_boundingbox=44,
                    resolution=0.5)
file.rename("./output","./output_no_currents_yes_depth")
```

## Generic data interpolation with no constraints - Land data interpolation

This function does not constraints for the interpolation. It uses the MCMC model to drive the final estimate around the prior interpolation. It could be used for quickly interpolating land data. The output is constituted of the same three files types described above.

Example:

```         
output<-bimac_onland(punctual_data_file,
                moving_average_points=1, 
                fast_solving=T,
                min_x_boundingbox=-0.25,
                max_x_boundingbox=12,
                min_y_boundingbox=34.75,
                max_y_boundingbox=44,
               resolution=0.5)
file.rename("./output","./output_no_currents_no_depth")
```

## Data preparation notes

The following steps are useful for retrieving and preparing current-velocity data:

1 - Download oceanic currents' components from Copernicus,e.g.:

<https://data.marine.copernicus.eu/product/GLOBAL_REANALYSIS_PHY_001_026/download?dataset=global-reanalysis-phy-001-026-grepv1-uv-monthly>

E.g., download a 1deg NetCDF file;

2 - Clip the file to the extent of the area to analyse, through QGIS;

3 - Save the components as two ASC files;

4 - Open the ASC files and check that the NO DATA value is -9999, otherwise run a string substitution of the NO DATA value to -9999;

5 - Download observation data from Argo: <https://dataselection.euro-argo.eu> Select delayed-mode only (if available) in the download panel to get higher quality data.
