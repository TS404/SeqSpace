#' Convert multiple sequence alignment to property matrix
#'
#' @param MSA        The reference MSA to be numericised
#' @param res.prop   The table of residue properties to use to numericsise the MSA
#'
#' @return Generates an object of class "Sequence Alignment Residue Properties" (SARP), providing the scales, numericised matrix generated from the MSA. Data is scaled within each property type. The details of the PAC output components are as follows:
#' \itemize{
#'  \item {MSA}             {Input MSA as a matrix of characters}
#'  \item {res.prop}        {Input residue property table as a data frame}
#'  \item {MSA.num.stack}   {Numericised MSA with each property as a separate list item}
#'  \item {MSA.num.wide}    {Numericised MSA with all properties in a single mastix}
#'  \item {MSA.scale.stack} {scaled, numericised MSA with each property as a separate list item}
#'  \item {MSA.scale.wide}  {scaled, numericised MSA with all properties in a single mastix}
#'  \item {prop.means}      {Means of each sequence property for the whole MSA}
#'  \item {prop.vars}       {Variances of each sequence property for the whole MSA}
#'  \item {seq.names}       {List of all sequence names}
#'  \item {seq.len}         {Length of sequences including gaps (i.e. length of whole MSA)}
#' } 
#' 
#' @export
#' @examples
#' data(example_MSA)
#' data(residue_properties)
#' SARP <- prepare_MSA(example_MSA,residue_properties)
#'
#' @note The `prepare_MSA' function converts an MSA into a matrix of scaled residue properties. Scaling is performed within each property type. The final matrix can be simplified by PCA
#'

#########
# Setup #
#########

library("tidyr")
library("seqinr")

prepare_MSA <- function(MSA,res.prop,
                        cys       = TRUE){
  
  # Load sequence MSA
  # if already a data frame, convert to matrix
  # if raw fasta file, use seqinr to convert
  if (is.data.frame(MSA)){
    MSA       <- as.matrix(t(toupper(as.matrix(MSA))))
    seq.names <- rownames(MSA)
    seq.len   <- ncol(MSA)
  }
  if (!is.matrix(MSA)){
    MSA       <- data.frame(seqinr::read.fasta(MSA,set.attributes=FALSE))
    MSA       <- as.matrix(t(toupper(as.matrix(MSA))))
    seq.names <- rownames(MSA)
    seq.len   <- ncol(MSA)
  }
  
  ######################
  # Residue properties #
  ######################
  
  # Load Residue properties
  if (!is.data.frame(res.prop)){
    res.prop <- read.csv(res.prop)
  }
  
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
  if (cys==1){
    res.prop[ncol(res.prop)+1]                    <- 0
    colnames(res.prop)[ncol(res.prop)]            <- "CYS"
    res.prop[grep("C", rownames(res.prop)),"CYS"] <- 1
  }
  # Append presence/absence columns
  res.prop[ncol(res.prop)+1]         <- c(rep(1,nrow(res.prop)-1),0)
  colnames(res.prop)[ncol(res.prop)] <- "NOTGAP"
  
  res.props <- colnames(res.prop)
  res.avail <- row.names(res.prop)
  
  
  ################
  # Numericising #
  ################

  # numericise MSA based on res.prop
  MSA.num.tall <- res.prop[t(MSA),]
  # name data types
  rownames(MSA.num.tall) <- NULL
  sequence     <- rep(x = seq.names, each  = seq.len)
  residue      <- rep(x = 1:seq.len, times = length(seq.names))
  MSA.num.tall <- cbind(sequence, residue, MSA.num.tall)
  
  # tidy data into list of matrices
  MSA.num.stack <- NULL
  for (x in 1:length(res.props)) {
    col.names <- paste(1:seq.len,
                       rep(res.props[x],seq.len),
                       sep = ".")
    MSA.num.stack[[res.props[x]]] <- matrix(MSA.num.tall[,x+2],
                                            ncol     = seq.len,
                                            byrow    = TRUE,
                                            dimnames = list(seq.names, 
                                                            col.names))
  }
  
  # Also reflow into single wide matrix
  MSA.num.wide <- MSA.num.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.num.wide <- cbind(MSA.num.wide, MSA.num.stack[[res.props[x]]])
  }
  
  
  ###############
  # Filter data #
  ###############
  
  # Assume no filtering of all-zero columns for now
  # MSA.filtered <- MSA.num[colMeans(MSA.num)!=0]
  # Generate list of remaining column names
  #MSA.names <- as.matrix(colnames(MSA.num.wide))
  #colnames(MSA.names) <- "name"
  #MSA.names <- tidyr::separate(data.frame(MSA.names),"name",
  #                         into   = c("resn","prop"),
  #                         sep    = "\\.",
  #                         remove = 1)
  #MSA.names <- t(MSA.names)
  
  
  ###############
  # Propscaling #
  ###############
  
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
  
  # replace gaps (currently "NA") with column average
  # create function
  na.colmean<-function(x){
    x[is.na(x)=="TRUE"] <- mean(as.matrix(x),na.rm = 1)
    x
  }
  # for each property of MSA.num.stack, apply function to each matrix comlumn
  for (x in 1:length(res.props)) {
    MSA.scale.stack[[x]] <- apply(MSA.scale.stack[[x]],2,na.colmean)
  }
  
  # Reflow into singe wide matrix for PCA
  MSA.scale.wide <- MSA.scale.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.scale.wide <- cbind(MSA.scale.wide, MSA.scale.stack[[x]])
  }
  
  
  #################
  # Ready for PCA #
  #################
  
  # Generate output object
  list(MSA             = MSA,
       res.prop        = res.prop,
       MSA.num.stack   = MSA.num.stack,
       MSA.num.wide    = MSA.num.wide,
       MSA.scale.stack = MSA.scale.stack,
       MSA.scale.wide  = MSA.scale.wide,
       prop.means      = prop.means,
       prop.vars       = prop.vars,
       seq.names       = seq.names,
       seq.len         = seq.len            )
}

  
  
