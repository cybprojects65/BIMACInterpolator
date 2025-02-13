cat("#DATA INTERPOLATION STARTED#\n")
rm(list=ls())
library(data.table)
library(raster)

t0=Sys.time()
set.seed(123) #for repeatability

#set the working directory to the current R file directory to preserve the paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#####SET INPUT AND AUX FILES
input_data<-"./input_data/temperature_argo_north_atlantic.csv"
output_data<-"./temperature_argo_north_atlantic_interpolated_without_currents.csv"
exclude_zeros<-F
longitude_column<-"longitude"
latitude_column<-"latitude"
value_column<-"temperature"
analysis_depth=-1
resolution=0.5 #degrees of the analysis resolution
#if the bimac distribution is multimodal perhaps the modelled quantity is not subject to advection-diffusion
#in this case, we suggest to take the distribution up to the maximum and set the other to -9999
#the following parameter can be used for such task
strict_interpolation<-F 

#####PROCESS START
#read the input file
data_to_interpolate<-as.data.frame(fread(input_data))
#select only the necessary columns
data_to_interpolate<-data_to_interpolate[,c(longitude_column,latitude_column,value_column)]
#recall the BIMAC functions
source("./BIMAC_functions.R")

#set the bounding box of the analysis
bounding_box_minX<-min(data_to_interpolate[longitude_column],na.rm = T)
bounding_box_maxX<-max(data_to_interpolate[longitude_column],na.rm = T)
bounding_box_minY<-min(data_to_interpolate[latitude_column],na.rm = T)
bounding_box_maxY<-max(data_to_interpolate[latitude_column],na.rm  = T)
bbox <- c(xmin = bounding_box_minX, xmax = bounding_box_maxX, ymin = bounding_box_minY, ymax = bounding_box_maxY)

#temporary aux file to be deleted at the end of the process
punctual_data_file<-"./tmp_bimac_data.csv"

#omit rows with at least one NA in one of the columns
data_to_interpolate<-na.omit(data_to_interpolate)

if (exclude_zeros==T){
  cat(" excluding zeros from the interpolation\n")
  data_to_interpolate<-data_to_interpolate[which(data_to_interpolate[value_column]>0),]
}else{
  cat(" including zeros in the interpolation\n")
}

write.csv(file=punctual_data_file,row.names = F,x = data_to_interpolate)

cat("#INTERPOLATING...\n")

output<-bimac_noadvection(punctual_data_file,
                          analysis_depth=analysis_depth,
                          moving_average_points=1, 
                          fast_solving=T,
                          min_x_boundingbox=bounding_box_minX,
                          max_x_boundingbox=bounding_box_maxX,
                          min_y_boundingbox=bounding_box_minY,
                          max_y_boundingbox=bounding_box_maxY,
                          resolution=resolution)

cat("#INTERPOLATION FINISHED...\n")
cat("THE OUTPUT HAS BEEN WRITTEN TO ./output/\n")
#retrieve the final bimac output
bimac_output<-"./output/BIMAC_interpolation.asc"
asc_file<-raster(bimac_output)
xseq<-seq(from=bounding_box_minX,to=bounding_box_maxX,by=resolution)
yseq<-seq(from=bounding_box_minY,to=bounding_box_maxY,by=resolution)
data_interpolated<-expand.grid(x = xseq, y = yseq) #combine the x and y coordinate to generate pairs
names(data_interpolated)<-c("x","y")
grid_values<-raster::extract(x=asc_file,y=data_interpolated,method='simple')
data_interpolated[,c(value_column)]<-grid_values

na_values_bimac<-which(is.na(data_interpolated[,c(value_column)]))
if (length(na_values_bimac)>0)
  data_interpolated<-data_interpolated[-na_values_bimac,]

par(mfrow = c(1, 2))
bimac_density<-density(data_interpolated[,c(value_column)])
plot(bimac_density$y,type='l')
points(bimac_density$y,col='red')
bimac_derivative<-bimac_density$y[2:length(bimac_density$y)]-bimac_density$y[1:(length(bimac_density$y)-1)]
sign_bimac_derivative<-as.numeric(sign(bimac_derivative))
first_descend<-min(which(sign_bimac_derivative==-1))-1
plot(bimac_density$y[1:first_descend],type='l',col='blue')
bimac_lower_limit<-bimac_density$x[first_descend]
cat("Maximum of the BIMAC distribution is",bimac_lower_limit,"\n")

if(strict_interpolation)
  data_interpolated[,c(value_column)][which(data_interpolated[,c(value_column)]>bimac_lower_limit)]<--9999

cat("Cleaning the cache..\n")
file.remove(punctual_data_file)
file.remove(" r2ssb.bug")
# Define the pattern for files: "<alphanumericcode>_gebco.csv"
pattern <- "^[A-Za-z0-9]+_gebco_30sec_8.asc$"
# List all files in the current directory
files <- list.files(path = "./", pattern = pattern, full.names = TRUE)
# If a matching file is found, delete it
if (length(files) > 0) {
  file.remove(files)
}


#save the output
cat("Saving the output to file..\n")
write.csv(file=output_data,x=data_interpolated,row.names = F)
t1=Sys.time()
cat("Process elapsed time:\n")
print(t1-t0)
cat("#DATA INTERPOLATION FINISHED#\n")
