new <- file.choose()
SAPCA <- PCAout

addsequence <- function(SAPCA,new){
  
  # Convert fasta 
  torotateSARP <- prepare_MSA(new)

  # Rotate the numericised data into the PCA space
  rotated <- scale(torotateSARP$num.scale.wide, SAPCA$PCA$center, SAPCA$PCA$scale) %*% SAPCA$PCA$rotation

  ###################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###################
  ###################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###################
  ################### Need what gaps were replaced with ###################
  ###################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###################
  ###################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###################

}









# Choose centres to rotate
# e.g. t(M$parameters$mean)->rotated
read.csv(file.choose()) -> torotate
# Convert first column to true 'row names'
names               <- c(as.matrix(torotate[,1]))  # convert first row of mydata to list via matrix
torotate            <- torotate[,2:ncol(torotate)] # remove original first column from 'mydata' 
row.names(torotate) <- names                       # name rows with 'names' list
rm (names)                                         # cleanup

# Rotate a row of data into the new PCA space
scale(torotate, pca$center, pca$scale) %*% pca$rotation -> rotated

# Rotate a row of data back out of PCA space
# NB: the data to be rotated must have the same numnber of dimensions as the rotation matrix
t(t(rotated %*% t(pca$rotation)) + pca$center)->reset

# Plot
plot3d(rotated,
    #   col    = c(k2$cluster,     # colour by clusters
    #              rep(5,1955)),   # ... and randoms
       #radius = 0.4,             # sphere radius if using spheres
       type   = "p",              # "p" is points, "s" is spheres
       size   = 4)                # point size   

# save
write.csv(reset, file=file.choose()) # sizes for spherical clusters






#####
# Stack randoms under real PCA data
rbind(pca2,rotated)->stacked
# Plot
plot3d(stacked,
          col    = c(k2$cluster,     # colour by clusters
                     rep(5,1955)),   # ... and randoms
       #radius = 0.4,             # sphere radius if using spheres
       type   = "p",              # "p" is points, "s" is spheres
       size   = 4)                # point size   