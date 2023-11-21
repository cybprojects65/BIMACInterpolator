rm(list=ls(all=TRUE))

source("BIMAC_functions.R")
currents_u_file = "./input_spatial_layers/oceanic_currents_u_wm.asc"
currents_v_file = "./input_spatial_layers/oceanic_currents_v_wm.asc"
punctual_data_file = "./observations/temperature_argo.csv"

output<-bimac_currents_depthboundaries(punctual_data_file,
                                       currents_u_file,
                                       currents_v_file,
                                       analysis_depth=-1,
                                       moving_average_points=1, 
                                       fast_solving=T, 
                                       sd_advection_equation=0.1
)

file.rename("./output","./output_yes_currents_yes_depth")

output<-bimac_nodepthboundaries(punctual_data_file,
                                currents_u_file,
                                currents_v_file, 
                                moving_average_points=1, 
                                fast_solving=T, 
                                sd_advection_equation=0.1
)

file.rename("./output","./output_yes_currents_no_depth")

output<-bimac_onland(punctual_data_file,
                     moving_average_points=1, 
                     fast_solving=T,
                     min_x_boundingbox=-0.25,
                     max_x_boundingbox=12,
                     min_y_boundingbox=34.75,
                     max_y_boundingbox=44,
                     resolution=0.5)

file.rename("./output","./output_no_currents_no_depth")

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
