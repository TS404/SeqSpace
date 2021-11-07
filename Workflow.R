#########
# Setup #
#########

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

# # Browse for labels file #OR# labels <- read.csv(file.choose())
# labels <- read.csv("C:\\Users\\T\\OneDrive\\0-Sequences\\2-PCA\\0-Raw data and scalers\\0 - Labels.csv")
# # Convert first column to true 'row names'
# names             <- c(as.matrix(labels[,1]))  # convert first row of mydata to list via matrix
# labels            <- labels[,2:ncol(labels)]   # remove original first column from 'mydata' 
# row.names(labels) <- names                     # name rows with 'names' list
# rm (names)                                     # cleanup


##############
# File names #
##############

# Numericise MSA
name      <- "VPE"
folder    <- "C:\\Users\\T\\OneDrive\\0-Sequences\\4-Non-defensins\\VPE"
res.prop1 <- "C:\\Users\\T\\OneDrive\\0-Sequences\\2-PCA\\0-Raw data and scalers\\Amino acid properties.csv"
# res.prop1 at https://github.com/TS404/SeqSpace/tree/master/data 

############
# Analysis #
############

setwd(folder)

SAPCA <- PCA_MSA (MSA = "VPE FULL LIST.fa", res.prop = res.prop1, cys=0, clusterPCs = 3, clusters = 1:20)
plot_3Dclusters(SAPCA, plotPCs = 1:3, labels = "Consensus")
plot_modelfit(SAPCA)
plot_loadings_matrix(SAPCA,PC=1,magnitude=FALSE)
plot_network(SAPCA, PC = 1:10, linkage = 20)

# SAPCA.C <- PCA_MSA (MSA = "Structureonly.fa", res.prop = res.prop1, cys=0, clusterPCs = 20, clusters = 1:8)
# plot_3Dclusters(SAPCA.C, plotPCs = 1:3, col = SAPCA$seq.space.clusters$classification,write = "test.3d.")
# plot_modelfit(SAPCA.C)
# 
# 
# SAPCA.F <- PCA_MSA (MSA = "FADSonly.fa", res.prop = res.prop1, cys=0, clusterPCs = 20, clusters = 1:8)
# plot_3Dclusters(SAPCA.F, plotPCs = 1:3, labels="gi303275396mine")
# plot_modelfit(SAPCA.F)

########
# Save #
########
newdir <- (paste(name,"PCAresults",sep="."))
dir.create(newdir)
setwd(paste(folder,newdir,sep = ""))
# Plots
plot_3Dclusters(SAPCA, labels = "Consensus",write = paste(name,"3Dplot.",sep="."))

svg(paste(name,"Scree.svg",sep="."))
  plot_scree(SAPCA)
dev.off()

svg(paste(name,"ModelFit.svg",sep="."))
  plot_modelfit(SAPCA)
  1
  0
dev.off()

i=1
svg(paste(name,".Loadings.PC",i,".svg",sep=""),width = 12, height = 1.5)
  plot_loadings_matrix(SAPCA,PC=i,topload=5,magnitude=FALSE)
dev.off()

i=2
svg(paste(name,".Loadings.PC",i,".svg",sep=""),width = 12, height = 1.5)
  plot_loadings_matrix(SAPCA,PC=i,topload=5,magnitude=FALSE)
dev.off()

i=3
svg(paste(name,".Loadings.PC",i,".svg",sep=""),width = 12, height = 1.5)
  plot_loadings_matrix(SAPCA,PC=i,topload=5,magnitude=FALSE)
dev.off()

i=4
svg(paste(name,".Loadings.PC",i,".svg",sep=""),width = 12, height = 1.5)
  plot_loadings_matrix(SAPCA,PC=i,topload=5,magnitude=FALSE)
dev.off()

i=5
svg(paste(name,".Loadings.PC",i,".svg",sep=""),width = 12, height = 1.5)
  plot_loadings_matrix(SAPCA,PC=i,topload=5,magnitude=FALSE)
dev.off()
 
# Data
saveRDS(SAPCA,file=paste(name,"PCA.RDS",sep="."))

# Export CSV files
for(i in 1:5){
   write.csv (topload(SAPCA,PC = i,n = 80), file = paste(name,".Loadings.PC",i,".csv",sep=""),row.names = FALSE)
}

write.csv (SAPCA$seq.space.PCA$loadings,           file = paste(name,"Loadings.All.csv",sep="."))
write.csv (SAPCA$seq.space.PCA$coordinates,        file = paste(name,"Coordinates.csv",sep=".")) 
clust.sum <- data.frame(SAPCA$seq.space.clusters$classification)
colnames(clust.sum) <- "cluster"
clust.sum[["certainty"]] <- apply(SAPCA$seq.space.clusters$likelihoods,1,max)
write.csv (clust.sum, file = paste(name,"Clusters.csv",sep="."))

setwd(folder)


###########################################
# Only structurally characterised columns #
###########################################

# Columns for Serpin alignment
notgap <- SAPCA$numerical.alignment$MSA["X134436.1.392",]!="-"


Structureregion  <- 185:694
Structurecolumns <- rep(FALSE,SAPCA$numerical.alignment$aln.len)
Structurecolumns[Structureregion] <- notgap[Structureregion]

tosave <- NULL
for (i in 1:3){
  temp <- loadingtable(SAPCA,
                       seq = SAPCA$numerical.alignment$MSA["X134436.1.392",],
                       magnitude = TRUE,
                       PC = i)
  tosave <- cbind(tosave,as.matrix(temp$annotated["sum",]))
}
 
write.csv(10*tosave[Structurecolumns,],paste(name,".Structureloadings.PC",PC,"(Pymol).csv", sep=""))




#######################
# Subcluster workflow #
#######################
subclusters=NULL
for (i in 1:SAPCA$clusters$G){
  subclusters[[i]] <- subPCA_MSA (SAPCA, cluster = i, clusterPCs = 5, clusters = 1:8)
}

sub.num <- 4
plot_3Dclusters(plotPCs = 1:3, subclusters[[sub.num]],
                col       = subset(labels$Label.simple.num,        subset=SAPCA$clusters$classification==sub.num), # colours list
                radius    = subset(0.075+labels$Antifungal.num/15, subset=SAPCA$clusters$classification==sub.num)) # sphere radius list)


# Generate graph
plot3d(subset(SAPCA$PCA$x,                               subset=SUB),
       col       = subset(labels$Label.simple.num,        subset=rownames(labels)==rownames(subcluster.1$PCA$x)), # colours list
       radius    = subset(0.075+labels$Antifungal.num/15, subset=rownames(labels)==rownames(subcluster.1$PCA$x)), # sphere radius list
       type      = "s",                                                # "p" is points, "s" is spheres
       specular  = "black")


#########
# Other #
#########

# Useful info from Mclust
SAPCA$cluster$parameters$mean             # centers of clusters
SAPCA$cluster$parameters$variance$scale   # size of clusters
SAPCA$cluster$modelName                   # optimal model

