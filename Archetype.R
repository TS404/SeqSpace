#########
# Setup #
#########

archetype <- function(SAPCA,
					  dims    = 3,
					  cluster = 1){
	# The number of dimensions that you care about the loadings for reconstruction
	# The cluster to use to create archetype
	if (is.numeric(cluster)){
		SUB <- SAPCA$seq.space.clusters$classification==cluster
	}else{
		SUB <- cluster
	}

	# The properties of the available residues
	# The PCA loadings for each residue property
	res.prop <- SAPCA$numerical.alignment$res.prop[1:(nrow(SAPCA$numerical.alignment$res.prop)-2),] # remove X and gap as options 
	arc.load <- SAPCA$seq.space.PCA$rotation[,1:dims]
	aln.len  <- SAPCA$numerical.alignment$aln.len

	# Gaps preference for subset
	consensus.gaps <- colMeans(subset(SAPCA$numerical.alignment$MSA.num.stack$NOTGAP, subset=SUB))


	######################
	# Residue properties #
	######################

	# Scale residue properties so that they can be searhed inside a unit circle around archetype properties similarly scaled
	res.prop.scale <- scale(res.prop)
	res.props      <- colnames(res.prop)
	res.avail      <- row.names(res.prop)


	######################
	# Archetype sequence #
	######################

	# Average residue properties of chosen subset
	arc.prop.num <- array(dim 	   = c(aln.len, 
	                                   length(res.props)),
	                      dimnames = list(NULL,
	                                      res.props))

	for (p in 1:length(res.props)) {
	  arc.prop.num[,p] <- colMeans(subset(SAPCA$numerical.alignment$MSA.num.stack[[p]],subset=SUB), na.rm=TRUE)
	}

	# is this the correct thing to be scaling? do I need to insert the typscaling step in the middle somewhere????????????????????????????
	# Scale sequence archetype requirements to same scale as the residue properties
	arc.prop.scale <- scale(arc.prop.num,
	                        center = attr(res.prop.scale,"scaled:center"),
	                        scale  = attr(res.prop.scale,"scaled:scale"))


	##########################
	# Loadings to prioritise #
	##########################

	# Positivised loadings for each property, summed for the number of PCA dimensions of interest
	arc.load.sum <- matrix(rowSums(sqrt(arc.load^2)),
						   nrow     = 194,
						   dimnames = list(1:aln.len,
										   res.props))


	############
	# Distance #
	############

	# Now we have the following:
	## properties of avaiable residues (res.prop.scale)
	## requred archetype sequence properties (arc.prop.scale)
	## the relative importance of each property to the PCA space (arc.load.sum)
	## overall distance is a^2 = b^2 + c^2 ... n^2
	# To choose the most apporpriate residue for each sequence archetype position we must find the distance of each position to each possible residue and identify the closest. PCA loading is used to search a spheroid, rather than spheerw 

	res.dist           <- array(dim=c(length(res.avail),
	                                   aln.len))
	rownames(res.dist) <- res.avail
	res.dist.load      <- res.dist # make empty copy

	# Distance is sqrt of all (arcproperty - resproperty)^2
	# Which residue is closest in a sphere around arcproperties
	for (x in 1:aln.len) {
	  temp.dis     <- t(arc.prop.scale[x,] - t(res.prop.scale))     # distance in each dimension
	  res.dist[,x] <- as.matrix(sqrt(rowSums(temp.dis^2,na.rm=1))) # root sum of squares for overall distance
	}

	# Loading-weighted distance is sqrt of all arcloading*(archpropery - resproperty)^2
	# Which residue is closest in a ovoid (determined by loadings) around arcproperties
	for (x in 1:aln.len) {
	  temp.dis          <- t(arc.prop.scale[x,] - t(res.prop.scale))        # distance in each dimension
	  temp.load         <- t(t(temp.dis) * sqrt(arc.load.sum[x,]))         # distances scaled by rooted loading weight
	  res.dist.load[,x] <- as.matrix(sqrt((rowSums(temp.load^2,na.rm=1)))) # root sum of squares
	}


	####################################
	# Select closest available residue #
	####################################

	# for each resuldue in res.dist, which residues has the lowest distance?
	# ignore those which should be gaps

	arc.seq.gapped <- NULL
	for (x in 1:aln.len) {
	  arc.seq.gapped[x] <- names(which.min(res.dist.load[,x]))
	}
	arc.seq.gapped[consensus.gaps<=0.5] <- "-"

	arc.seq <- paste(grep("[^-]",(arc.seq.gapped),value=TRUE),collapse="")

	# record the distances of the closest available residues
	arc.seq.dist      <- apply(res.dist,2,min)
	arc.seq.dist.load <- apply(res.dist.load,2,min)

	arc.seq.dist     [consensus.gaps<=0.5] <- NA
	arc.seq.dist.load[consensus.gaps<=0.5] <- NA

	#############
	# Consensus #
	#############

	# Consensus sequence for comparison
	con.seq.gapped <- seqinr::consensus(subset(SAPCA$numerical.alignment$MSA, subset=SUB))
	con.seq        <- paste(grep("[^-]",(con.seq.gapped),value=TRUE),collapse="")

	##########
	# Output #
	##########

	# Generate output object
	list(arc.seq          = arc.seq,
		 arc.seq.gapped   = arc.seq.gapped,
		 arc.prop.num     = arc.prop.num,
		 arc.load         = arc.load,
		 con.seq          = con.seq,
		 con.seq.gapped   = con.seq.gapped,
		 call = list(SAPCA   = SAPCA,
		 			 dims    = dims,
		 			 cluster = cluster)
		)
}


# Also, maybe arc.seq.dist(.load) could be scaled 0-1?