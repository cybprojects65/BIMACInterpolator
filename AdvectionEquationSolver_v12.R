rm(list=ls(all=TRUE))
library(raster)
library ( R2jags )
library ( coda )
library(plyr)
library(dplyr)

#ASC files definitions
currents_u_file = "./input_spatial_layers/oceanic_currents_u_wm.asc"
currents_v_file = "./input_spatial_layers/oceanic_currents_v_wm.asc"
depth_file = "./input_spatial_layers/gebco_30sec_8_clip_025deg_wm.asc"

#currents_u_file = "./input_spatial_layers/oceanic_currents_u_1deg.asc"
#currents_v_file = "./input_spatial_layers/oceanic_currents_v_1deg.asc"
#depth_file = "./input_spatial_layers/gebco_30sec_8_clipped_res_1.asc"

punctual_data_file = "./observations/temperature_argo.csv"

#OUTPUT files will be: 
prior_data_output = gsub(pattern = ".csv",replacement = "_advsolver_IW_prior.asc",x = punctual_data_file)
posterior_data_output = gsub(pattern = ".csv",replacement = "_advsolver_posterior.asc",x = punctual_data_file)

#Function to transform a raster file into a matrix
asc_to_matrix<-function(asc_raster_file){

  asc_file<-raster(asc_raster_file)
  min_x_in_raster<-asc_file@extent[1]
  max_x_in_raster<-asc_file@extent[2]
  min_y_in_raster<-asc_file@extent[3]
  max_y_in_raster<-asc_file@extent[4]
  resolution<-res(asc_file)[1]
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
  
  return(values)
}

#Function to check the consistency between the velocity and depth files
check_file_consistency<-function(currents_u_file,currents_v_file,depth_file){
  u<-raster(currents_u_file)
  v<-raster(currents_v_file)
  d<-raster(depth_file)
  r_u<-resolution<-res(u)[1]
  r_v<-resolution<-res(v)[1]
  r_d<-resolution<-res(d)[1]
  if ((r_u!=r_v) || (r_u!=r_d)){
    cat("ERROR: misaligned resolutions in files\n")
    return(F)
  }
  else{
    if ((u@extent[1]!=v@extent[1]) || (u@extent[1]!=d@extent[1])){
      cat("ERROR: misaligned in left longitude between files\n")
      return(F)
    }
    if ((u@extent[2]!=v@extent[2]) || (u@extent[2]!=d@extent[2])){
      cat("ERROR: misaligned in right longitude between files\n")
      return(F)
    }
    if ((u@extent[3]!=v@extent[3]) || (u@extent[3]!=d@extent[3])){
      cat("ERROR: misaligned in lower latitude between files\n")
      return(F)
    }
    if ((u@extent[4]!=v@extent[4]) || (u@extent[4]!=d@extent[4])){
      cat("ERROR: misaligned in upper latitude between files\n")
      return(F)
    }
    return(T)
  }
  
}

#Step 1 - Create matrices out of the input files
cat("Reading data files\n")
currents_u_matrix<-asc_to_matrix(currents_u_file)
currents_v_matrix<-asc_to_matrix(currents_v_file)
depth_matrix<-asc_to_matrix(depth_file)

if (!check_file_consistency(currents_u_file,currents_v_file,depth_file)){
  stop("ERROR: HETEROGENEOUS INPUT FILES - PLEASE PROVIDE U-V-Depth SPATIALLY ALIGNED FILES AT THE SAME RESOLUTION\n")
}else{
  cat("   Input files are consistent\n")
}

#Step 2 - Retrieve points on land as positive bathymetry points
cat("Retrieving land\n")
land<-which(as.vector(t(depth_matrix))>=0)
if (length(land)==0) {
  land<-c()
  cat("NOTE: NO LAND LOCATIONS ARE PRESENT\n")
}

#Step 3 - Retrieve observation points, filter those falling in the velocity bounding box, and assign them to a regular grid filled with NA values
cat("Retrieving observation points\n")
asc_raster_file<-currents_u_file
asc_file<-raster(asc_raster_file)
min_x_in_raster<-asc_file@extent[1]
max_x_in_raster<-asc_file@extent[2]
min_y_in_raster<-asc_file@extent[3]
max_y_in_raster<-asc_file@extent[4]
resolution<-res(asc_file)[1]
#create the grid
xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
grid_of_points<-expand.grid(x = xseq, y = yseq) #combine the x and y coordinate to generate pairs
grid_val<-grid_of_points
grid_val$value<-0
#read the observation points
punctual_data<-read.csv(file = punctual_data_file)
names(punctual_data)<-c("x","y","value")
ypoints<-unique(grid_of_points$y)
xpoints<-unique(grid_of_points$x)
ncol_r<-length(xpoints)
nrow_r<-length(ypoints)
#delete points falling out of the velocity bounding box
cat("Filtering observation points\n")
punctual_data_filt<-punctual_data[which(punctual_data$x>min_x_in_raster & punctual_data$x<max_x_in_raster),]
punctual_data_filt<-punctual_data_filt[which(punctual_data_filt$y>min_y_in_raster & punctual_data_filt$y<max_y_in_raster),]
if (dim(punctual_data_filt)[1]==0){
  stop("ERROR: NO POINT PRESENT IN THE ANALYSIS AREA\n")
}
#calculate the matrix indices of the observations
punctual_data_filt$x_i<-round((punctual_data_filt$x-min_x_in_raster)/resolution)+1
punctual_data_filt$y_i<-round((max_y_in_raster-punctual_data_filt$y)/resolution)+1
n_punctual_data<-dim(punctual_data_filt)[1]
data_matrix<-matrix(nrow = nrow_r,ncol = ncol_r,data = NA)
#assign the observations to matrix elements
cat("Transforming observations into a matrix\n")
for (k in 1:nrow(punctual_data_filt)){
  row_p<-punctual_data_filt [k,]
  data_matrix[row_p$y_i,row_p$x_i]<-row_p$value
}
min_real_value<-min(data_matrix,na.rm = T)
max_real_value<-max(data_matrix,na.rm = T)

#Step 4 - calculate the average proximity between the points to assess a radius for prior value averaging
cat("Calculating average proximity between the points\n")
min_distances<-sapply(1:n_punctual_data, function(i){
  x<-punctual_data_filt$x[i]
  y<-punctual_data_filt$y[i]
  dist<-sqrt( ( (x-punctual_data_filt$x[-i])*(x-punctual_data_filt$x[-i]) ) + ( (y-punctual_data_filt$y[-i])*(y-punctual_data_filt$y[-i]) ) ) 
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
      index_in_vector<-((i-1)*ncol(data_matrix))+j
      #assign -9999 to land points
      if ( (length(land)>0) && (index_in_vector%in%land) )
        avg<- -9999
      else{
        #calculate the inverse distance with respect to all points
        dist<-sqrt( ( (x-punctual_data_filt$x)*(x-punctual_data_filt$x) ) + ( (y-punctual_data_filt$y)*(y-punctual_data_filt$y) ) ) 
        inv_dist<-1/(dist)
        #select the distances within the proximity range
        good_points<-which(dist<=upper_limit_min_distances)
        if (length(good_points)>0){
          #if there are points within the proximity range, calculate the inverse-distance weighted average and assign it to the current cell
          avg<-sum((inv_dist[good_points]*punctual_data_filt$value[good_points]))/sum(inv_dist[good_points])
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
  for (i in 1:nrow_r){
    for (j in 1:ncol_r){
      element<-data_matrix_prefilled[i,j]
      #if the current element is a hole then process it
      if (is.na(element)){
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
          #avg<-mean(submatrix[valid_indices]) simple average - obsolete
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
  #increase the number of iterations
  iterations<-iterations+1
}
cat("\n")

#Step 7 - save the inverse weighted map to a file
cat("Saving the inverse weighted map\n")
ro <- raster(ncol=ncol_r, nrow=nrow_r)
length(values(ro))
data_matrix_prefilled_vector<-as.vector(t(data_matrix_prefilled))
data_matrix_prefilled_vector_to_save<-data_matrix_prefilled_vector
if (length(land)>0)
  data_matrix_prefilled_vector_to_save[land]<-NA
values(ro)<-data_matrix_prefilled_vector_to_save
extent(ro)<-extent(asc_file)
NAvalue(ro)<- -9999
writeRaster(ro, filename=prior_data_output, format="ascii",overwrite=TRUE)

#Step 8 - solve the advection equation and recalculate the grid-matrix values
cat("Solving advection equation...\n")
#data_matrix_postfilled is the posterior values' matrix
data_matrix_postfilled<-matrix(nrow = nrow_r,ncol = ncol_r,data = 1)
data_matrix_postfilled<-as.vector(data_matrix_postfilled)
N<-length(data_matrix_postfilled)
#calculate the standard deviation of the observation values
prior_sd<-sd(data_matrix,na.rm=T)
data_vector<-as.vector(t(data_matrix))
if (length(which(!is.na(data_vector)))==0){
  stop("ERROR IN THE DATA: NO REAL OBSERVATION PRESENT IN THE AREA!")
}
if (length(which(is.na(data_vector)))==0){
  stop("ERROR IN THE DATA: THERE ARE NO POINTS TO RECONSTRUCT!")
}
#record the indexes of the observation values in the data matrix
is_true_value<-!is.na(data_vector)
real_values_idx<-which(!is.na(data_vector))
#record the indexes of the non-observation values to recalculate
values_to_estimate_idx<-which(is.na(data_vector))

#define a zero-filled array as the value of the advection function
zeros<-integer(length(is_true_value))
#define velocity vectors to use in the advection function
ux<-as.vector(t(currents_u_matrix))
uy<-as.vector(t(currents_v_matrix))

#define the values for which calculating advection is possible
valid_index_for_advection<-matrix(nrow = nrow_r,ncol = ncol_r,data = 1)
# exclude the the first colum
valid_index_for_advection[,1]<-0
# exclude the first row
valid_index_for_advection[1,]<-0
# exclude the the last colum
valid_index_for_advection[,ncol_r]<-0
# exclude the the last row
valid_index_for_advection[nrow_r,]<-0
#transform the matrix into an array for easier modelling
valid_index_for_advection_vector<-as.vector(t(valid_index_for_advection))

#exclude land points from the data
land_no_values<-c()
if (length(land)>0){
  #exclude land points from the non-observation values indexes
  values_to_estimate_idx<-values_to_estimate_idx[-which(values_to_estimate_idx%in%land)]
  #set land fluid velocity to 0
  ux[land]<-0
  uy[land]<-0
  #exclude land points from being involved in advection calculation, otherwise coastal underestimation can occur
  valid_index_for_advection_vector[land]<-0
  #exclude possible observation values falling on land from the land indexes - these should be used as likelihood data
  true_in_land<-which(land%in%real_values_idx)
  if (length(true_in_land)>0){
    land_no_values<-land[-true_in_land]
  }else
    land_no_values<-land
  #set all land values to zero
  data_matrix_prefilled_vector[land_no_values]<-0  
}
#select valid indexes on which calculate the advection equation 
valid_index_for_advection_vector_idx<-which(valid_index_for_advection_vector!=0)

#Set the input parameters of the model
#prior_sd = standard deviation of the observation values
#data_matrix_prefilled_vector = priors' vector, zero for land values
#min_real_value = minimum observation value
#max_real_value = maximum observation value
#ux = horizontal current velocity
#uy = vertical current velocity
#zeros = all-zero results of the advection equation
#ncol_r = number of columns in the data matrix used to select advection terms correctly
#resolution = the data spatial resolution
#real_values_idx = vector indices of the real observations
#values_to_estimate_idx = vector indices of the values to estimate
#valid_index_for_advection_vector_idx = indices of values valid for the advection equation
#land_no_values = indices of land locations excluding those containing real observations
jags.data <- list ("data_matrix_prefilled_vector",
                   "min_real_value","max_real_value",
                   "ux","uy",
                   "prior_sd",
                   "zeros",
                   "ncol_r",
                   "resolution",
                   "real_values_idx",
                   "values_to_estimate_idx",
                   "valid_index_for_advection_vector_idx",
                   "land_no_values")
#P = posterior distribution values
jags.params <- c("P")

#Note: we use differential equations like https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/advection.pdf
Model = "
model {
#priors

#force posterior distribution values to 0
for (g in land_no_values){
  P[g]<-0
}

#initialise the individual-location posterior distribution values in real observation on to the real-values' range
for (g in real_values_idx){
  P[g] ~ dunif (
      (min_real_value),
      (max_real_value))
}

#fit the individual distribution values to real observations in the real-value locations - this accounts for error in the observations
for (g in real_values_idx){
  invsigma[g] <- pow(0.1,-2)
  data_matrix_prefilled_vector[g] ~ dnorm(P[g],invsigma[g])
}

#inizialise the other individual-location posterior distribution values to inverse-weighted estimated observations
for (g in values_to_estimate_idx){
  invsigma_na[g] <- pow(prior_sd,-2)#pow(prior_sd,-2)
  P[g] ~ dnorm (data_matrix_prefilled_vector[g],invsigma_na[g])
}

#relate the posterior distribution through the advection equation
for (k in valid_index_for_advection_vector_idx){
  #if the terms involve land points set them to 0
  ux_dpsi_dx[k]<-ifelse(P[k+1]==0 || P[k-1]==0, 0, 
    ux[k]*(P[k+1]-P[k-1])/(2*resolution)
  )
  uy_dpsi_dy[k]<-ifelse(P[k-ncol_r]==0 || P[k+ncol_r]==0, 0,
    uy[k]*(P[k-ncol_r]-P[k+ncol_r])/(2*resolution)
    )
  #advection term: it is set to 0 if at least one term is 0
  advection[k]<-ifelse(ux_dpsi_dx[k]==0 || uy_dpsi_dy[k]==0, 
                0, 
                ux_dpsi_dx[k]+uy_dpsi_dy[k])
  invsigma_adv[k] <- pow(0.1,-2)
  #fit the equation to 0
  zeros[k] ~ dnorm(advection[k],invsigma_adv[k])
}

}"

#write the BUGS model to a file
JAGSFILE =" r2ssb.bug "
cat (Model , file = JAGSFILE )
Nchains = 1 # number of Markov chains - to account for non - ergodic convergence
Nburnin = 100 # burn -in iterations - n. of initial iterations to discard
Niter = 1000 # total n. of iterations
Nthin = 10 # thinning - take every 10 samples to lower the dependency among the samples
#Run the Gibbs sampling
jagsfit <- jags ( data = jags.data , working.directory=NULL , inits =NULL , jags.params,
                  model.file = JAGSFILE , n.chains=Nchains, n.thin=Nthin , n.iter=Niter , n.burnin=Nburnin )
#retrieve the posterior vector
dmp<-jagsfit$BUGSoutput$sims.list$P
#the posterior vector has one column for each P element, samples from the gibbs sampling are contained in each column
#the mean extract the optimal value according to the Monte Carlo integration theorem
dmpC<-colMeans(dmp)
#set land points to NA
if(length(land)>0)
  dmpC[land]<-NA


#Step 9 - distribution smoothing
cat("Smoothing\n")
data_matrix_postfilled<-matrix(dmpC, nrow = nrow_r,ncol = ncol_r,byrow = T)
smooth = T
if (smooth){
  ma<-1
  data_matrix_postfilled_smoothed<-sapply(1:nrow_r, function(i){
    #take ma-elements around each point in a row
    min_i<-max(0,(i-ma))
    max_i<-min(nrow_r,(i+ma))
    smoothed<-integer(ncol_r) #inizialise with zeros
    smoothed<-sapply(1:ncol_r, function(j){
      #take ma-elements around each point along the column
      min_j<-max(1,(j-ma))
      max_j<-min(ncol_r,(j+ma))
      #extract the submatrix
      sub_matrix<-data_matrix_postfilled[min_i:max_i,min_j:max_j]
      #transform the submatrix into a vector and delete NA values (land points)
      sub_vector<-as.vector(sub_matrix)
      valid_indices_subvector<-which(!is.na(sub_vector))
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
  #transpose since sapply returns rows as columns
  data_matrix_postfilled_smoothed<-t(data_matrix_postfilled_smoothed)
  #substitute the smoothed matrix to the posterior matrix
  data_matrix_postfilled<-data_matrix_postfilled_smoothed
  
  cat("Saving smoothed map...\n")
}

#Step 10 - save the posterior probability matrix
ro_p <- raster(ncol=ncol_r, nrow=nrow_r)
values_vec_post<-as.vector(t(data_matrix_postfilled))
if (length(land)>0)
  values_vec_post[land]<--9999
values(ro_p)<-values_vec_post
extent(ro_p)<-extent(asc_file)
NAvalue(ro_p)<- -9999
writeRaster(ro_p, filename=posterior_data_output, format="ascii",overwrite=TRUE)
prior<-as.vector(t(data_matrix_prefilled))
posterior<-as.vector(t(data_matrix_postfilled))

#evaluate the difference between posterior and prior probabilities
if (length((land)>0)){
  prior[land]<-NA
  posterior[land]<-NA
}
diff<-posterior-prior
if (length((land)>0)){
  diff<-abs(diff[-which(is.na(posterior-prior))])
  diff_rel<-diff/abs(posterior[-land])
}else{
  diff<-abs(diff)
  diff_rel<-diff/abs(posterior)
}
cat("max discrepancy",max(diff),"\n")
cat("mean relative discrepancy",100*mean(diff_rel),"%\n")
cat("sd discrepancy",sd(diff),"\n")
cat("Done.\n")  
