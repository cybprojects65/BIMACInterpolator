#BIMAC interpolation algorithms - created by Gianpaolo Coro - gianpaolo.coro@cnr.it
rm(list=ls(all=TRUE))
start_time <- Sys.time()
library(raster)
library ( R2jags )
library ( coda )
library(plyr)
library(dplyr)
library(digest)

#DEFAULT EXAMPLE: Global scale temperature
#ASC files definitions
punctual_data_file = "./observations/temperature_argo.csv"

#MODEL PARAMETERS
depth_file = "gebco_30sec_8.asc"
analysis_depth_min<--12000 #discard all points over this limit
analysis_depth_max<--1
automatic_extent<-T
resolution<-1
moving_average_points<-1
smooth = T

#OUTPUT files will be: 
prior_data_output = "./output/BIMAC_IDW_prior_no_advection.asc"

if (moving_average_points==0)
  smooth = F

#Step 1 - Create matrices out of the input files
cat("Reading data file\n")
punctual_data<-read.csv(file = punctual_data_file)
names(punctual_data)<-c("x","y","value")


if (automatic_extent){
  min_x_in_raster<-min(punctual_data$x,na.rm = T)
  max_x_in_raster<-max(punctual_data$x,na.rm = T)
  min_y_in_raster<-min(punctual_data$y,na.rm = T)
  max_y_in_raster<-max(punctual_data$y,na.rm = T)
}

adapt_depth<-function(depth_file, resolution, min_x_in_raster, min_y_in_raster, max_x_in_raster, max_y_in_raster){
  caching_string<-paste0(resolution,";",min_x_in_raster,";",min_y_in_raster,";",max_x_in_raster,";",max_y_in_raster)
  sha<-sha1(caching_string)
  cachedfile<-paste0(sha,"_",depth_file)
if (!file.exists(cachedfile)){
    
  asc_file<-raster(depth_file)
  xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
  yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
  grid_of_points<-expand.grid(x = xseq, y = yseq) #combine the x and y coordinate to generate pairs
  grid_values<-extract(x=asc_file,y=grid_of_points,method='simple') #extract raster values for the observations and the grid
  grid_val<-grid_of_points
  grid_val$value<-grid_values
  ypoints<-unique(grid_of_points$y)
  xpoints<-unique(grid_of_points$x)
  ncol_r<-length(xpoints)
  nrow_r<-length(ypoints)
  values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
  row_counter<-1
  for (y_c in 1:(nrow_r)){
    yp<-ypoints[y_c]
    row_rast<-grid_val[which(grid_val$y == yp),]
    row_rast<-row_rast[order(row_rast$x),]
    values[(nrow_r-row_counter+1),]<-row_rast$value[1:(ncol_r)]
    row_counter<-row_counter+1
  }
  save(values,file=cachedfile)
}else{
  cat("...Loading cached file",cachedfile,"\n")
  load(cachedfile)
}
  return(values)
}

depth_matrix<-adapt_depth(depth_file, resolution, min_x_in_raster, min_y_in_raster, max_x_in_raster, max_y_in_raster)

#Step 2 - Retrieve points on land as positive bathymetry points
cat("Retrieving land\n")
land<-which( (as.vector(t(depth_matrix))>analysis_depth_max) | (as.vector(t(depth_matrix))<analysis_depth_min))

if (length(land)==0) {
  land<-c()
  cat("NOTE: NO LAND LOCATIONS ARE PRESENT\n")
}else if (length(land)==length(as.vector(t(depth_matrix)))) {
  stop("ERROR: No point at the given depth is present in the area, please change the analysis depth input") 
}

#Step 3 - building data grid
cat("Building observation grid\n")
#create the grid
xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
grid_of_points<-expand.grid(x = xseq, y = yseq) #combine the x and y coordinate to generate pairs
grid_val<-grid_of_points
grid_val$value<-0
#read the observation points
ypoints<-unique(grid_of_points$y)
xpoints<-unique(grid_of_points$x)
ncol_r<-length(xpoints)
nrow_r<-length(ypoints)

#delete points falling out of the velocity bounding box
cat("Filtering observation points\n")
punctual_data_within_area<-punctual_data[which(punctual_data$x>min_x_in_raster & punctual_data$x<max_x_in_raster),]
punctual_data_within_area<-punctual_data_within_area[which(punctual_data_within_area$y>min_y_in_raster & punctual_data_within_area$y<max_y_in_raster),]
if (dim(punctual_data_within_area)[1]==0){
  stop("ERROR: NO POINT PRESENT IN THE ANALYSIS AREA\n")
}
#calculate the matrix indices of the observations
punctual_data_within_area$x_i<-round((punctual_data_within_area$x-min_x_in_raster)/resolution)+1
punctual_data_within_area$y_i<-round((max_y_in_raster-punctual_data_within_area$y)/resolution)+1
n_punctual_data<-dim(punctual_data_within_area)[1]
data_matrix<-matrix(nrow = nrow_r,ncol = ncol_r,data = NA)
#assign the observations to matrix elements
distinct_observation_pairs<-distinct(data.frame(x_i=punctual_data_within_area$x_i,y_i=punctual_data_within_area$y_i))
distinct_observation_pairs_xyv<-list()

cat("Transforming observations into a matrix through IDW\n")
inserted_observations<-1
for (pair in 1:nrow(distinct_observation_pairs)){
      xi = distinct_observation_pairs[pair,]$x_i
      yi = distinct_observation_pairs[pair,]$y_i
      pk<-punctual_data_within_area[which(punctual_data_within_area$x_i==xi & punctual_data_within_area$y_i == yi),]
      dist<-sqrt( ( (xi-pk$x)*(xi-pk$x) ) + ( (yi-pk$y)*(yi-pk$y) ) ) 
      if (length(which(dist==0))>0)
        dist[which(dist==0)] = resolution/0.1
      
      inv_dist<-1/(dist)
      avg<-sum((inv_dist*pk$value)/sum(inv_dist))
      if (yi<=nrow(data_matrix) && xi<=ncol(data_matrix)){
        data_matrix[yi,xi]<-avg
        #index to coordinates
        x_m<-((xi-1)*resolution)+min_x_in_raster
        y_m<-max_y_in_raster-((yi-1)*resolution)
        distinct_observation_pairs_xyv[[inserted_observations]]<-data.frame(x_m=x_m,y_m=y_m,value=avg)
        inserted_observations<-inserted_observations+1
      }
}
min_real_value<-min(data_matrix,na.rm = T)
max_real_value<-max(data_matrix,na.rm = T)
distinct_observation_pairs_xyv <- ldply(distinct_observation_pairs_xyv, data.frame)

#Step 4 - calculate the average proximity between the points to assess a radius for prior value averaging
cat("Calculating average proximity between the points\n")
min_distances<-sapply(1:nrow(distinct_observation_pairs_xyv), function(i){
  x<-distinct_observation_pairs_xyv$x_m[i]
  y<-distinct_observation_pairs_xyv$y_m[i]
  ddx<-(x-distinct_observation_pairs_xyv$x_m[-i])
  ddy<-(y-distinct_observation_pairs_xyv$y_m[-i])
  dist<-sqrt( (ddx*ddx) + (ddy*ddy))
  return(min(dist))
},simplify = T)

#use the geometric mean
log_mean_min_dist<-mean(log(min_distances[which(min_distances>0)]))
log_sd_min_dist<-sd(log(min_distances[which(min_distances>0)]))
upper_limit_min_distances<-exp(log_mean_min_dist+1.96*log_sd_min_dist)
upper_limit_min_distances_index<-round(upper_limit_min_distances/resolution)
cat("Proximity range is",upper_limit_min_distances,"deg","(=",upper_limit_min_distances_index,"indices)","\n")

#Step 5 - calculate prior values of the quantity through inverse weighting
cat("Inverse weighting\n")
#data_matrix_prefilled is the matrix of prior values
data_matrix_prefilled<-data_matrix

on_land<-function(i,j,resolution,min_x_in_raster,data_matrix,land){
  x<-((j-1)*resolution)+min_x_in_raster
  y<-max_y_in_raster-((i-1)*resolution)
  index_in_vector<-((i-1)*ncol(data_matrix))+j
  if ( (length(land)>0) && (index_in_vector%in%land) )
    return(T)
  else
    return(F)
}

cat("   Filling gaps around the observations\n")
#for each NA element (=without observation) retrieve observation points within the proximity range; calculate the inverse-weighted values and assign the weighted average
data_matrix_prefilled<-sapply(1:nrow_r, function(i){
  data_matrix_prefilled_row<-vector()
  for (j in 1:ncol_r){
    element<-data_matrix[i,j]
    #if the element is not an observation then process it
    if (is.na(element)){
      #retrieve the x,y coordinates of the cell
      x<-((j-1)*resolution)+min_x_in_raster
      y<-max_y_in_raster-((i-1)*resolution)
      #assign -9999 to land points
      if (on_land(i,j,resolution,min_x_in_raster,data_matrix,land)){
        avg<- -9999
      }
      else{
        #calculate the inverse distance with respect to all points
        ddx<-(x-distinct_observation_pairs_xyv$x_m)
        ddy<-(y-distinct_observation_pairs_xyv$y_m)
        dist<-sqrt( (ddx*ddx) + (ddy*ddy))
        inv_dist<-1/(dist)
        #select the distances within the proximity range
        good_points<-which(dist<=upper_limit_min_distances)
        if (length(good_points)>0){
          #if there are points within the proximity range, calculate the inverse-distance weighted average and assign it to the current cell
          avg<-sum((inv_dist[good_points]*distinct_observation_pairs_xyv$value[good_points]))/sum(inv_dist[good_points])
        }else{ 
          #if no proximity point was found force a hole in the matrix
          avg<- NA
        }
      }
      data_matrix_prefilled_row[j]<-avg
    }else{
      data_matrix_prefilled_row[j]<-element
    }
  }
  return(data_matrix_prefilled_row)
},simplify = T)
#adjust the matrix after the sapply, which returns the rows as columns
data_matrix_prefilled<-t(data_matrix_prefilled)

#Step 6 - Fill the holes' matrix holes through an iterative gap-filling process
cat("   Averaging gaps\n")
cat("Iterating.. ")
still_elements_to_fill<-T
iterations<-1
#iterate until no holes are present
while (still_elements_to_fill){
  still_elements_to_fill<-F
  elements_to_fill<-0
  elements_to_fill_prev<-0
  for (i in 1:nrow_r){
    for (j in 1:ncol_r){
      element<-data_matrix_prefilled[i,j]
      #if the current element is a hole then process it
      if (is.na(element) && !(on_land(i,j,resolution,min_x_in_raster,data_matrix,land))){
        avg<-NA
        #retrieve all cells within a proximity distance (this time taken as a matrix index)
        di<-upper_limit_min_distances_index#1 #distance as an index
        i_0 = max(1,(i-di))
        i_1 = min(nrow_r,(i+di))
        j_0 = max(1,(j-di))
        j_1 = min(ncol_r,(j+di))
        i_c = ((i_1-i_0)/2)+1
        j_c = ((j_1-j_0)/2)+1
        
        submatrix<-as.vector(data_matrix_prefilled[i_0:i_1,j_0:j_1])
        submatrix_dist<-as.vector(data_matrix_prefilled[i_0:i_1,j_0:j_1])
        #calculate the inverse-distance values depending on how far the points are from the current point
        submatrix_vector_index<-1
        for (gy in j_0:j_1){
          dy = abs(gy-j)*resolution
          for (gx in i_0:i_1){
            dx = abs(gx-i)*resolution
            id = 1/((dx*dx)+(dy*dy))
            if (is.na(data_matrix_prefilled[gx,gy]) || (data_matrix_prefilled[gx,gy]==-9999) ){
                id = 1
            }
            submatrix_dist[submatrix_vector_index]<-id
            submatrix_vector_index=submatrix_vector_index+1
          }
        }
        
        #retrieve only non-hole and non-land points
        valid_indices<-which(!is.na(submatrix) & submatrix!=-9999)
        if (length(valid_indices)>0){
          #if there are points around fill the hole with the inverse-distance weighted average
          sm_valid<-submatrix[valid_indices]
          sm_valid_d<-submatrix_dist[valid_indices]
          avg<-sum(sm_valid_d*sm_valid)/sum(sm_valid_d)
        }else{
          #if there hole remains then alert the loop that holes are still present
          still_elements_to_fill<-T
          #increase the number of holes
          elements_to_fill<-elements_to_fill+1
        }
        data_matrix_prefilled[i,j]<-avg
      }
    }
  }
  cat(iterations,"(",elements_to_fill,"); ")
  if ((elements_to_fill_prev>0) && (elements_to_fill<elements_to_fill_prev)){
    elements_to_fill_prev<-elements_to_fill
  }else{
    upper_limit_min_distances_index<-upper_limit_min_distances_index+1
  }
  #increase the number of iterations
  iterations<-iterations+1
}
cat("\n")

if (smooth){
#smooth
  data_matrix_prefilled_smoothed<-sapply(1:nrow_r, function(i){
    #take moving_average_points-elements around each point in a row
    min_i<-max(1,(i-moving_average_points))
    max_i<-min(nrow_r,(i+moving_average_points))
    smoothed<-integer(ncol_r) #inizialise with zeros
    smoothed<-sapply(1:ncol_r, function(j){
      #take moving_average_points-elements around each point along the column
      min_j<-max(1,(j-moving_average_points))
      max_j<-min(ncol_r,(j+moving_average_points))
      #extract the submatrix
      sub_matrix<-data_matrix_prefilled[min_i:max_i,min_j:max_j]
      #transform the submatrix into a vector and delete NA values (land points)
      sub_vector<-as.vector(sub_matrix)
      valid_indices_subvector<-which(!is.na(sub_vector) & sub_vector!=-9999)
      #if points exist report the submatrix mean for the current matrix element; otherwise report -9999
      if (length(valid_indices_subvector)>0){
        avg<-mean(sub_matrix[valid_indices_subvector])
      }
      else
        avg<--9999
      return(avg)
    },simplify = T)
    
    return(smoothed)
  },simplify = T)
data_matrix_prefilled<-t(data_matrix_prefilled_smoothed)
}
cat("Saving the inverse weighted map\n")
ro <- raster(ncol=ncol_r, nrow=nrow_r)
length(values(ro))
data_matrix_prefilled_vector<-as.vector(t(data_matrix_prefilled))
data_matrix_prefilled_vector_to_save<-data_matrix_prefilled_vector

if (length(land)>0){
  data_matrix_prefilled_vector_to_save[land]<--9999
}
values(ro)<-data_matrix_prefilled_vector_to_save
extent(ro)<-extent( min_x_in_raster+(resolution/2),
                    max_x_in_raster-(resolution/2),
                    min_y_in_raster+(resolution/2),
                    max_y_in_raster-(resolution/2))
NAvalue(ro)<- -9999
writeRaster(ro, filename=prior_data_output, format="ascii",overwrite=TRUE)

end_time <- Sys.time()
cat("Done.\n")
cat("Elapsed.\n")
print(end_time-start_time)
