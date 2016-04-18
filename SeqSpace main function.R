#' Perform PCA on sequence property matrix
#'
#' @param SARP        The Sequence Alignment Residue Properties object produced by `prepare_PCA'
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

#library("rgl")
#library("rglwidget")
#library("ape")
#library("ggplot2")
#library("tidyr")
#library("seqinr")


PCA_MSA <- function(MSA,
                    res.prop,
                    cys        = TRUE,
                    clusterPCs = 5,
                    clusters   = 1:8){
  
  # for some reason, seem to need to explicitly load this library
  library(quietly = TRUE,"mclust")
  
  ###############
  # Prepare MSA #
  ###############

  # Load sequence MSA
  # if a matrix, can be used straight away
  # if raw fasta file, use seqinr to convert to data frame
  if (!is.matrix(MSA)){
    MSA       <- data.frame(seqinr::read.fasta(MSA,set.attributes=FALSE))
  }
  # if a data frame, convert to matrix
  if (is.data.frame(MSA)){
    MSA       <- as.matrix(t(toupper(as.matrix(MSA))))
  }

  # currrently deal with this issue at later step by replacing NA-only columns with 0s
      ## remove empty columns (all gaps) since they screw up later steps
      #if(length(which(colMeans(MSA=="-")==1))>=1){
      #  emptycolumn <- which(colMeans(MSA=="-")==1)
      #  MSA <- MSA[,-emptycolumn]
      #  print(paste("MSA column ",emptycolumn," was empty and was removed"))
      #}
  
  seq.names <- rownames(MSA)
  aln.len   <- ncol(MSA)  

  ##############################
  # Prepare residue properties #
  ##############################
  
  # Load residue properties
  # if a data frame, can be used straight away
  # if a matrix, convert to data frame 
  if (is.matrix(res.prop)){
    res.prop <- data.frame(res.prop)
  }
  # if a csv file, convert to data frame
  if (!is.data.frame(res.prop)){
    res.prop <- read.csv(res.prop)
  }
  
  # Format rownames and order
  row.names(res.prop) <- c(as.matrix(res.prop[,1]))           # first column to row names
  res.prop            <- res.prop[,-1]                        # remove original first column
  res.prop            <- res.prop[,order(colnames(res.prop))] # order columns alphabetically
  # Append X and gap rows (if X not already defined)
  if (length(grep("X", rownames(res.prop)))==0){
    res.prop          <- rbind(res.prop, colMeans(res.prop))
    rownames(res.prop)[nrow(res.prop)] <- "X"
  }
  res.prop            <- rbind(res.prop, rep(NA,ncol(res.prop)))
  rownames(res.prop)[nrow(res.prop)]   <- "-"
  # Add cystene columns (gaps = 0)
  if (cys==TRUE){
    res.prop[ncol(res.prop)+1]                    <- 0
    colnames(res.prop)[ncol(res.prop)]            <- "CYS"
    res.prop[grep("C", rownames(res.prop)),"CYS"] <- 1
  }
  # Append presence/absence columns
  res.prop[ncol(res.prop)+1]         <- c(rep(1,nrow(res.prop)-1),0)
  colnames(res.prop)[ncol(res.prop)] <- "NOTGAP"
  
  res.props <- colnames(res.prop)
  res.avail <- row.names(res.prop)
  
  
  ##################
  # Numericise MSA #
  ##################

  # Numericise MSA based on res.prop
  MSA.num.tall <- res.prop[t(MSA),]
  # Name data types
  rownames(MSA.num.tall) <- NULL
  sequence     <- rep(x = seq.names, each  = aln.len)
  residue      <- rep(x = 1:aln.len, times = length(seq.names))
  MSA.num.tall <- cbind(sequence, residue, MSA.num.tall)
  # Stack data into list of matrices
  MSA.num.stack <- NULL
  for (x in 1:length(res.props)) {
    col.names <- paste(1:aln.len,
                       rep(res.props[x],aln.len),
                       sep = ".")
    MSA.num.stack[[res.props[x]]] <- matrix(MSA.num.tall[,x+2],
                                            ncol     = aln.len,
                                            byrow    = TRUE,
                                            dimnames = list(seq.names, 
                                                            col.names))
  }
  # Also reflow into single wide matrix
  MSA.num.wide <- MSA.num.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.num.wide <- cbind(MSA.num.wide, MSA.num.stack[[res.props[x]]])
  }
  
  
  ############################
  # Scaling by property type #
  ############################
  
  # Take means and variances of each property type
  prop.means <- NULL
  prop.vars  <- NULL
  for (x in 1:length(res.props)) {
    prop.means[x] <- mean(MSA.num.stack[[x]],na.rm=1)
    prop.vars[x]  <- var(tidyr::gather(data.frame(MSA.num.stack[[x]]))[2],na.rm=1)
  }
  names(prop.means) <- res.props
  names(prop.vars)  <- res.props
  
  # Scale numericised MSA to prop.means and prop.vars
  MSA.scale.stack <- NULL
  for (x in 1:length(res.props)) {
    MSA.scale.stack[[res.props[x]]] <- (MSA.num.stack[[res.props[x]]]- prop.means[x]) /
                                        sqrt(prop.vars[x])
  }
  
  # Replace gaps (currently "NA") with column average
  # Create na.colmean function
  na.colmean<-function(x){
    x[is.na(x)=="TRUE"] <- mean(as.matrix(x),na.rm = 1)
    x
  }
  # For each property of MSA.num.stack, apply na.colmean function to each matrix comlumn
  for (x in 1:length(res.props)) {
    MSA.scale.stack[[x]] <- apply(MSA.scale.stack[[x]],2,na.colmean)
  }

  # Also reflow into singe wide matrix for PCA
  MSA.scale.wide <- MSA.scale.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.scale.wide <- cbind(MSA.scale.wide, MSA.scale.stack[[x]])
  }
    

  ##################
  # Alignment list #
  ##################

  numerical.alignment <- list(MSA             = MSA,
                              res.prop        = res.prop,
                              MSA.num.stack   = MSA.num.stack,
                              MSA.num.wide    = MSA.num.wide,
                              MSA.scale.stack = MSA.scale.stack,
                              MSA.scale.wide  = MSA.scale.wide,
                              prop.means      = prop.means,
                              prop.vars       = prop.vars,
                              seq.names       = seq.names,
                              aln.len         = aln.len)

  

  ##############################################################
  # Principal Component analysis of numericised sequence space #
  ##############################################################

  # Prepare data for PCA 
  toPCA <- numerical.alignment$MSA.scale.wide
  # If any columns contained gaps only, replace all valuses with 0
  toPCA[,is.na(colMeans(toPCA))]<-0
  
  if (cys==FALSE){
    toPCA <- toPCA[,grep("CYS", colnames(toPCA), invert = TRUE)]
  }
  
  # PCA of data with no extra scaling
  PCA.raw <- stats::prcomp(toPCA)

  seq.space.PCA <- list(coordinates = PCA.raw$x,
                        loadings    = PCA.raw$rotation,
                        centre      = PCA.raw$center,
                        scale       = PCA.raw$scale,
                        stdev       = PCA.raw$sdev)
            
  ################################
  # Finding model-based clusters #
  ################################
  # Mclust calculating multiple cluster numbsers and chooses most likely
  # Different models assume different (non-circular) cluster shapes
  # Different numbers of clusters
  # Manual - http://www.stat.washington.edu/research/reports/2012/tr597.pdf
  # Models - http://finzi.psych.upenn.edu/R/library/mclust/html/mclustModelNames.html
  
  # Using Mclust on PCA data
  clusters.raw <- mclust::Mclust(seq.space.PCA$coordinates[,1:clusterPCs], # which PCs to use to find clusters
                                 prior = priorControl(),         # starting values for BIC iterations
                                 G     = clusters)               # number of possible clusters to assess
  
  clusters.min                <- clusters.raw
  clusters.min$classification <- NULL
  clusters.min$G              <- NULL
  clusters.min$z              <- NULL
  clusters.min$call           <- NULL

  seq.space.clusters <- list(classification   = clusters.raw$classification,
                             optimal          = clusters.raw$G,
                             checked          = clusters,
                             likelihoods      = clusters.raw$z,
                             other            = clusters.min)


  ##########
  # Output #
  ##########

  list(numerical.alignment = numerical.alignment,
       seq.space.PCA       = seq.space.PCA,
       seq.space.clusters  = seq.space.clusters, 
       call                = list(MSA        = MSA,
                                  res.prop   = res.prop,
                                  cys        = cys,
                                  clusterPCs = clusterPCs,
                                  clusters   = clusters))
}







############
# Analysis #
############

topload <- function(SAPCA,
                    PC = 1,
                    n  = 20){
  
  names     <- do.call(rbind, strsplit(gsub("\\.",":",rownames(SAPCA$seq.space.PCA$loadings)), ':'))
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


loadtable <- function(SAPCA,
                      PC){

    data <- matrix(data     = sqrt(SAPCA$seq.space.PCA$loadings[,PC]^2),
                 nrow     = length(colnames(SAPCA$numerical.alignment$res.prop)),
                 byrow    = TRUE,
                 dimnames = list(colnames(SAPCA$numerical.alignment$res.prop),
                                   1:SAPCA$numerical.alignment$aln.len))
  sum  <- colSums(data)
  data <- rbind(SAPCA$numerical.alignment$MSA["gi303275396mine",], data, sum)  
}








#########
# PLOTS #
#########
plot_modelfit <- function(SAPCA){
  # Model plots     
  plot(SAPCA$seq.space.clusters$other$BIC) # Using plot on Mclust data reqires numbers to be added afterwards (not sure why)
}

plot_scree <- function(SAPCA){
  # Scree plot of component significance
  barplot(SAPCA$seq.space.PCA$stdev[1:15],               # first 15 principal components
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
    colour <- SAPCA$seq.space.clusters$classification
  }else{
    colour <- col
  }
  # Calculate radius size
  if (all(radius == "auto")){
    rad <- (range(SAPCA$seq.space.PCA$coordinates[,1])[2]-range(SAPCA$seq.space.PCA$coordinates[,1])[1])/100
  }else{
    rad <- radius
  }
  # Plot model-based clusters in 3D
  rgl::plot3d(SAPCA$seq.space.PCA$coordinates[,plotPCs],
              col      = colour,      # colour by clusters
              specular = "black",     # matte lighting
              type     = "s",         # "p" is points, "s" is spheres
              radius   = rad,         # sphere radius if using spheres
              size     = 4)           # point size
 
  for (NAME in labels){
    SUB = row.names(SAPCA$seq.space.PCA$coordinates)==NAME      # Label based on its row.name
    rgl::text3d(subset(SAPCA$seq.space.PCA$coordinates,subset=SUB), 
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


plot_loadings <- function (SAPCA){

  barnames <- 1:SAPCA$numerical.alignment$aln.len
  data     <- matrix(data  = sqrt(SAPCA$seq.space.PCA$loadings[,1]^2),
                     nrow  = length(colnames(SAPCA$numerical.alignment$res.prop)),
                     byrow = FALSE)
  
  barplot(data      = data,
          names.arg = barnames)
  
}