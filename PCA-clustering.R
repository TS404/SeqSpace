#' Perform PCA on sequence property matrix
#'
#' @param SARP        The Sequence Alignment Residue Properties object produced by `prepare_PCA'
#' @param labels      CSV file of labels of any data to be associated with sequences
#' @param clusterPCs  The number of principal components to be used to identify clusters of sequences (default 5)
#' @param clusters    The number of sequence clusters to search for (default = 1:8)
#'
#' @return Generates an object of class "Sequence Alignment Residue Properties" (SARP), providing the scales, numericised matrix generated from the MSA. Data is scaled within each property type. The details of the PAC output components are as follows:
#' \itemize{
#'  \item {MSA}             {Input MSA as a matrix of characters}
#'  \item {res.prop}        {Input residue property table as a data frame}
#'  \item {MSA.num.stack}   {Numericised MSA with each property as a separate list item}
#' } 
#' 
#' @export
#' @examples
#' data(example_MSA)
#' data(residue_properties)
#' SARP <- prepare_MSA(example_MSA,residue_properties)
#' PCA  <- PCA_MSA(SARP)
#'
#' @note The `PCA_MSA' function performs multidimensional scaling of a numericised MSA into a reduced space by PCA, to simplify it and aid interpretation.
#'

#########
# Setup #
#########

library("rgl")
library("rglwidget")
library("ape")
library("mclust")
library("ggplot2")


PCA_MSA <- function(SARP,
                    labels     = NULL,
                    clusterPCs = 5,
                    clusters   = 1:8){
  
  ###############
  # Import data #
  ###############
  
  if (!is.null(SARP$MSA.scale.wide)) {
    data <- SARP$MSA.scale.wide
  }
  if (is.matrix(SARP)) {
    data <- data
  }
  if (all(is.null(SARP$MSA.scale.wide),!is.matrix(SARP))) {
    data            <- read.csv(SARP)
    row.names(data) <- c(as.matrix(data[,1])) # name rows with 'names' list
    data            <- data[,-1]              # remove original first column from 'mydata'
  }
  
  if (occupancy == FALSE){
    data <- data[,-grep("NOTGAP", colnames(data))]
  }
  
  if (cys == FALSE){
    data <- data[,-grep("CYS", colnames(data))]
  }
  
  #######
  # PCA #
  #######
  # PCA of data with no extra scaling
  pca <- prcomp(data)
  
  ################################
  # Finding model-based clusters #
  ################################
  # Mclust calculating multiple cluster numbsers and chooses most likely
  # Different models assume different (non-circular) cluster shapes
  # Different numbers of clusters
  # Manual - http://www.stat.washington.edu/research/reports/2012/tr597.pdf
  # Models - http://finzi.psych.upenn.edu/R/library/mclust/html/mclustModelNames.html
  
  # Using Mclust on PCA data
  M <- mclust::Mclust(pca$x[,1:clusterPCs],         # which PCs to use to find clusters
                      prior = priorControl(),       # starting values for BIC iterations
                      G = clusters)                 # number of possible clusters to assess
  M$call<-NULL
  
  ##########
  # Output #
  ##########
  
  list(PCA      = pca,
       clusters = M, 
       SARP     = SARP)
}



############
# Analysis #
############

topload <- function(SAPCA,
                    PC = 1,
                    n  = 20){
  
  names     <- do.call(rbind, strsplit(gsub("\\.",":",rownames(SAPCA$PCA$rotation)), ':'))
  consensus <- seqinr::consensus(SAPCA$SARP$MSA)
  combined  <- cbind(names,consensus,SAPCA$PCA$rotation[,PC])
  sorted    <- combined[order(sqrt(SAPCA$PCA$rotation[,PC]^2),decreasing = TRUE),]
  rownames(sorted) <- NULL
  colnames(sorted) <- c("resn","property","consensus",paste("PC",PC,"_load",sep=""))
  if (n=="all") {
    sorted
  }else{
    head(sorted,n)
  }
}



#########
# PLOTS #
#########
plot_modelfit <- function(SAPCA){
  # Model plots     
  plot(SAPCA$clusters) # Using plot on Mclust data reqires numbers to be added afterwards (not sure why)
  1  # Baysian Information Criterion = fit of each model with different cluster numbers
  #2  # Coloured cluters
  #3  # Uncertainty (darker dots less certain)
  0  # exit
}

plot_scree <- function(SAPCA){
  # Scree plot of component significance
  barplot(SAPCA$PCA$sdev[1:15],               # first 15 principal components
          xlab = "Pricipal component",        # x label
          ylab = "Variance",                  # y label
          main = "Principal components")      # title
}

plot_3Dclusters <- function(SAPCA,
                            plotPCs = 1:3,
                            col     = "cluster",
                            radius  = "auto",
                            labels  = NULL,
                            write   = FALSE){
  if (all(col=="cluster")){
    colour <- SAPCA$clusters$classification
  }else{
    colour <- col
  }
  # Calculate radius size
  if (all(radius == "auto")){
    rad <- (range(SAPCA$PCA$x[,1])[2]-range(SAPCA$PCA$x[,1])[1])/100
  }else{
    rad <- radius
  }
  # Plot model-based clusters in 3D
  rgl::plot3d(SAPCA$PCA$x[,plotPCs],
              col      = colour,      # colour by clusters
              specular = "black",     # matte lighting
              type     = "s",         # "p" is points, "s" is spheres
              radius   = rad,         # sphere radius if using spheres
              size     = 4)           # point size
 
  for (NAME in labels){
    SUB = row.names(SAPCA$PCA$x)==NAME      # Label based on its row.name
    text3d(subset(SAPCA$PCA$x,subset=SUB), 
           text      = paste('---',NAME),   # data label text
           font      = 2,                   # bold
           color     = "black",             # colour
           adj       = -0.3)                # offset
  }
  
  # Write html for interactive data
  if (write!=FALSE){
    rglwidget::.writeWebGL(write)                          
  }
}


