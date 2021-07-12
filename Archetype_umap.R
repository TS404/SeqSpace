# Also, maybe arc.seq.dist(.load) could be scaled 0-1?
archetype.umap <- function(umap = fas.msa.umap,
                           hdbscan = fas.msa.hdbscan,
                           cluster = 1){
  # The number of dimensions that you care about the loadings for reconstruction
  # The cluster to use to create archetype
  if (is.numeric(cluster)){
    row.SUB <- hdbscan$cluster==cluster
  }else{
    row.SUB <- cluster
  }
  
  col.SUB <- colMeans(fas.msa.umap$data[row.SUB,grep(pattern = "NOTGAP", x = colnames(fas.msa.umap$data))])>=0
  differences <- (colMeans(fas.msa.umap$data[row.SUB,col.SUB]) - colMeans(fas.msa.umap$data[,col.SUB])) 
  hist(differences,breaks = 80)
}


