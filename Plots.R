############
# Analysis #
############

topload <- function(SAPCA,
                    PC = 1,
                    n  = 20){
  
  names     <- do.call(rbind, strsplit(gsub("\\.",":",rownames(SAPCA$PCA$rotation)), ':'))
  consensus <- seqinr::consensus(SAPCA$numerical.alignment$MSA)
  combined  <- cbind(names,consensus,SAPCA$seq.space.PCA$loadings[,PC])
  sorted    <- combined[order(sqrt(SAPCA$seq.space.PCA$loadings[,PC]^2),decreasing = TRUE),]
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
  # Model plots of either clusters or subclusters
  if(!isnull(SAPCA$seq.space.clusters$others$BIC)){
    data <- SAPCA$seq.space.clusters$others$BIC
  }
  if(!isnull(SAPCA$seq.space.subclusters$others$BIC)){
    data <- SAPCA$seq.space.subclusters$others$BIC
  }
  plot(data)
}

plot_scree <- function(SAPCA){
  # Scree plot of component significance of either clusters or subclusters
  if(!isnull(SAPCA$seq.space.PCA$stdev)){
    data <- SAPCA$seq.space.PCA$stdev
  }
  if(!isnull(SAPCA$seq.space.cluster.PCA$stdev)){
    data <- SAPCA$seq.space.cluster.PCA$stdev
  }

  barplot(data[1:15],                    # first 15 principal components
          xlab = "Pricipal component",   # x label
          ylab = "Variance",             # y label
          main = "Principal components") # title
}

plot_3Dclusters <- function(SAPCA,
                            plotPCs = 1:3,
                            col     = "cluster",
                            radius  = "auto",
                            labels  = NULL,
                            write   = FALSE){
  # Plot 3D scatterplot of either clusters or subclusters
  if(!isnull(SAPCA$seq.space.PCA$coordinates)){
    data <- SAPCA$seq.space.PCA$coordinates
  }
  if(!isnull(SAPCA$seq.space.cluster.PCA$coordinates)){
    data <- SAPCA$seq.space.cluster.PCA$coordinates
  }

  # Datapoint colours
  if (all(col=="cluster")){
    colour <- SAPCA$seq.space.clusters$classification
  }else{
    colour <- col
  }
  # Datapoint sizes
  if (all(radius == "auto")){
    rad <- (range(SAPCA$seq.space.PCA$coordinates[,1])[2]-range(SAPCA$seq.space.PCA$coordinates[,1])[1])/100
  }else{
    rad <- radius
  }

  # Plot model-based clusters in 3D
  rgl::plot3d(data[,plotPCs],
              col      = colour,      # colour by clusters
              specular = "black",     # matte lighting
              type     = "s",         # "p" is points, "s" is spheres
              radius   = rad,         # sphere radius if using spheres
              size     = 4)           # point size
 
  for (NAME in labels){
    SUB = row.names(SAPCA$seq.space.PCA$coordinates)==NAME      # Label based on its row.name
    text3d(subset(SAPCA$seq.space.PCA$coordinates,subset=SUB), 
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