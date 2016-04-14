#############################
# Gather necessary packages #
#############################

# Activate packages
library(quiet=TRUE,"rgl")
library(quiet=TRUE,"ape")
library(quiet=TRUE,"cluster")
library(quiet=TRUE,"mclust")
library(quiet=TRUE,"gplots")
library(quiet=TRUE,"ggplot2")
# Define colour scheme
colours<-palette(c("blue",       #1
                   "red",        #2
                   "yellow",     #3
                   "white",      #4
                   "darkgrey",   #5
                   "maroon",     #6
                   "orange",     #7
                   "black",      #8
                   "purple",     #9
                   "darkgreen")) #10

subPCA_MSA <- function(SAPCA,
                       cluster    = 1,
                       clusterPCs = 5,
                       clusters   = 1:8){
  
  #################
  # Define subset #
  #################
  
  # Define subset for graph
  SUB = SAPCA$clusters$classification==cluster
        #labels$M.AA.type.pc10.8Gs==1|2|3

  ###########################
  # SUBSET PCA and clusters #
  ###########################
  
  # Create subset of data and PCA with no extra scaling
  subPCA  <- prcomp(data.frame(subset(SAPCA$SARP$MSA.scale.wide,SUB)))                       
  
  ################################
  # Finding model-based clusters #
  ################################
  # Mclust calculating multiple cluster numbsers and chooses most likely
  # Different models assume different (non-circular) cluster shapes
  # Different numbers of clusters
  # Manual - http://www.stat.washington.edu/research/reports/2012/tr597.pdf
  # Models - http://finzi.psych.upenn.edu/R/library/mclust/html/mclustModelNames.html
  
  # Using Mclust on PCA data
  M <- mclust::Mclust(subPCA$x[,1:clusterPCs],         # which PCs to use to find clusters
                      prior = priorControl(),       # starting values for BIC iterations
                      G = clusters)                 # number of possible clusters to assess
  M$call<-NULL
  
  ##########
  # Output #
  ##########
  
  list(PCA      = subPCA,
       clusters = M, 
       SAPCA    = SAPCA)
}
