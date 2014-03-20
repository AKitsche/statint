#Function to display interaction contrast matrices 
plotcontrast <- function(cmat){
  if(!is.matrix(cmat) || !is.numeric(cmat)){stop("cmat must be a matrix with numeric entries")}
  if(is.null(rownames(cmat))){rownames(cmat)<-paste("c", 1:nrow(cmat), sep="")}
  
  dcm <- melt(cmat, varnames = c("comparison", "treatment"), value.name="coefficient" )
  dcm$comparison <- factor(dcm$comparison, levels=rev(rownames(cmat)))
  pp <- ggplot(dcm, aes(y=comparison, x=treatment)) + 
    theme_grey() +
    geom_tile(aes(fill=coefficient), colour="black") +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) 
  
  return(pp)
}
