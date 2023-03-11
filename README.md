# BIMAC - Bayesian Interpolation Model with Advection-diffusion Constraint
An advection equation solver for oceanic parameter spatial interpolation based on Jags and Gibbs sampling.

The software is an all-in-one script (BIMAC.R). 
The essential input parameters to change are:

currents_u_file = "./input_spatial_layers/oceanic_currents_u_1deg.asc"

currents_v_file = "./input_spatial_layers/oceanic_currents_v_1deg.asc"

punctual_data_file = "./observations/temperature_argo.csv"

analysis_depth<--1

![](https://github.com/cybprojects65/JagsOceanicSpatialInterpolator/blob/main/global_scale_example.png)
*Figure 1. Example of global scale temperature layer estimated from observation data from the Argo network.* 

## Notes on data preparation

1 - Download oceanic currents' components from Copernicus,e.g.:

https://data.marine.copernicus.eu/product/GLOBAL_REANALYSIS_PHY_001_026/download?dataset=global-reanalysis-phy-001-026-grepv1-uv-monthly

2 - Download a 1deg NetCDF file

3 - Clip the file to the area extent

4 - Save the components as two ASC files

5 - Open the ASC files and check that the NO DATA value is -9999, otherwise run a string substitution of the NO DATa value to -9999

6 - To compare the results with DIVA, clip the NetCDF extent to the same extent of the velocity files

7 - Export the the netcdf bands to ASC files

8 - Change xll and yll ASC values to the velocity files values plus half resolution (e.g., -180.5 -71.5 if the velocity files were -180.0, -71.0). This correct a known bug of the DIVA files.

9 - Download observation data from Argo: https://dataselection.euro-argo.eu

10 - Select delayed-mode only (if available) in the download panel.
