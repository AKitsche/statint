#function iacontrast() to create named direct Kronecker interaction contrasts

# fa, fb: two factor variables, each with at least two levels
# typea, typeb: character string, one of those as described for type in contrMat (multcomp), ignored if cma, cmb are given, respectively. If provided, typea and typeb are ignored.
# cma, cmb: contrast matrices, with as many columns as there are levels in fa, fb, respectively (before or after dropping unused levels, see droplevels below).
# abbrevnames=NULL: a named list of arguments to be passed to abbreviate to abbreviate the factor level names for the construction of names for the interaction contrasts
# sep = " - ", character string to be used to separate different levels of the cellmeans in comparisons in the contrasts names
# cw = "," character string, to be used when collapsing original factor levels within a factor (pooling contrasts)
# cb = "," character string, to be used when collapsing original factor levels to cellmean levels in the contrasts names
# droplevels=TRUE: logical, indicating whether unused levels of the input factors fa, fb should be dropped. If TRUE, the number of columns in cma, cmb, should fit the number of used levels. If FALSE, the number of columns in cma, cmb shoudl fit the number of level in the leve?l attribute of the factors.

# POSSIBLE ADD. ARGUMNETS (NOT YET IMPLEMENTED):
# orderby=NULL: CURRENTLY NOT USED a,b, whether the output contrasts should be additionally ordered acc. to factor a, b 
# basea, baseb: CURRENTLY NOT USED, to be passed to contrMat, if typea, typeb="Dunnett", respectively.


# VALUE: a list with two elements:
# fab is the cellmeans factor (with one level for each level combination of fa, fb)
# cmab the interaction contrast matrix, with columns and naming of contrast fitting the level order in fab.


iacontrast <- function(fa, fb, typea="Dunnett", typeb="Dunnett", cma=NULL, cmb=NULL, droplevels=TRUE, abbrevnames=NULL, orderby=NULL, sep=" - ", cw=",", cb=":", method="at")
{
  method <- match.arg(method, choices=c("at", "cellmeans"))
  
  if(!is.factor(fa)){warning("Argument fa has been coerced to a factor!"); FA <- as.factor(fa)}else{if(droplevels){FA<-droplevels(fa)}else{FA<-fa}}
  if(!is.factor(fb)){warning("Argument fb has been coerced to a factor!"); FB <- as.factor(fb)}else{if(droplevels){FB<-droplevels(fb)}else{FB<-fb}}
  
  nga <- table(FA)
  ngb <- table(FB)
  
  if(length(nga)<2){if(droplevels){stop("Factor specified in fa should have at least two levels (after dropping unused levels).")}else{stop("Factor specified in fa should have at least two levels")}}
  if(length(ngb)<2){if(droplevels){stop("Factor specified in fb should have at least two levels (after dropping unused levels).")}else{stop("Factor specified in fb should have at least two levels")}}
  
  
  if(is.null(cma)){CMA<-contrMat(n=nga, type=typea)}else{
    if(!is.matrix(cma)){stop("Argument cma is not a matrix.")}
    if(!is.numeric(cma)){stop("Argument cma must be a matrix with numeric entries.")}
    if(ncol(cma)!=length(nga)){if(droplevels){stop("Argument cma must have as many columns as there are levels in the factor specified in fa, after dropping unused levels.")}else{stop("Argument cma must have as many columns as there are levels in the factor specified in fa.")}}
    rsa <- abs(rowSums(cma))
    if(any(rsa > 10 * .Machine$double.eps)){warning("At least one contrast (row) in cma does not sum to zero.")}
    CMA <- cma
  }
  
  
  if(is.null(cmb)){CMB<-contrMat(n=ngb, type=typeb)}else{
    if(!is.matrix(cmb)){stop("Argument cmb is not a matrix.")}
    if(!is.numeric(cmb)){stop("Argument cmb must be a matrix with numeric entries.")}
    if(ncol(cmb)!=length(ngb)){if(droplevels){stop("Argument cmb must have as many columns as there are levels in the factor specified in fb, after dropping unused levels.")}else{stop("Argument cmb must have as many columns as there are levels in the factor specified in fb.")}}
    rsb <- abs(rowSums(cmb))
    if(any(rsb > 10 * .Machine$double.eps)){warning("At least one contrast (row) in cmb does not sum to zero.")}
    CMB <- cmb
  }
  
  nca<-nrow(CMA)
  ncb<-nrow(CMB)
  
  
  if(!is.null(abbrevnames))
  {
    
    if(!is.list(abbrevnames)){stop("Argument abbrevnames must be a list.")}
    
    LABBA<-abbrevnames
    LABBA$names.arg <- levels(FA) 
    NAMESA<-do.call("abbreviate", args=LABBA)
    
    LABBB<-abbrevnames
    LABBB$names.arg <- levels(FB) 
    NAMESB<-do.call("abbreviate", args=LABBB)
    
  }else{
    NAMESA<-levels(FA)
    NAMESB<-levels(FB)
  }
  
  
  scma <- t(apply(CMA, MARGIN=1, sign))
  scmb <- t(apply(CMB, MARGIN=1, sign))
  
  all0a <- apply(scma, MARGIN=1, function(x){all(x==0)})
  all0b <- apply(scmb, MARGIN=1, function(x){all(x==0)})
  
  if(any(all0a)){stop("Argument cma contains at least one contrast (row) containing only zeros.")}
  if(any(all0b)){stop("Argument cmb contains at least one contrast (row) containing only zeros.")}
  
  
  posnamesa <- apply(scma, MARGIN=1, function(x){paste(NAMESA[x==1], collapse=cw)})
  posnamesb <- apply(scmb, MARGIN=1, function(x){paste(NAMESB[x==1], collapse=cw)})
  
  negnamesa <- apply(scma, MARGIN=1, function(x){paste(NAMESA[x==-1], collapse=cw)})
  negnamesb <- apply(scmb, MARGIN=1, function(x){paste(NAMESB[x==-1], collapse=cw)})
  
  
  FAB <- FA:FB
  
  if(any(table(FAB)==0)){warning("Factors specified in fa and fb are not completely cross classified. For some level combinations in fab there are no observations.")}
  
  CMAB <- kronecker(CMA, CMB) 
  
  posnamesav <- rep(posnamesa, each=ncb)
  negnamesav <- rep(negnamesa, each=ncb)
  
  posnamesbv <- rep(posnamesb, times=nca)
  negnamesbv <- rep(negnamesb, times=nca)
  
  # If both matrices are ordinary contrast matrices:
  
  if(all(negnamesa!="")& all(negnamesb!=""))
  {
    switch(method,
           "at"={
             COMPNAMESA <- paste("(",paste(posnamesav, negnamesav, sep=sep), ")", sep="")
             COMPNAMESAB <- paste(paste("(",paste(COMPNAMESA, posnamesb, sep=cb), ")", sep=""), paste("(",paste(COMPNAMESA, negnamesb, sep=cb), ")", sep=""), sep=sep)
           },
           "cellmeans"={
             PNAMESAB <- paste( "(", paste( paste(posnamesav, posnamesbv, sep=cb) , paste(negnamesav, posnamesbv, sep=cb), sep=sep), ")", sep="")
             NNAMESAB <- paste( "(", paste( paste(posnamesav, negnamesbv, sep=cb) , paste(negnamesav, negnamesbv, sep=cb), sep=sep), ")", sep="")
             COMPNAMESAB <- paste(PNAMESAB, NNAMESAB, sep=sep)
           })
  }else{
    
    # If cma is not an ordinary contrast matrix (i.e., there are NO negative entries in cma)
    
    if(all(negnamesa=="")){
      
      switch(method,
             "at"={
               COMPNAMESB <- paste("(",paste(posnamesbv, negnamesbv, sep=sep), ")", sep="")
               COMPNAMESAB <- paste(COMPNAMESB, posnamesav, sep=cb)
             },
             "cellmeans"={
               PNAMESAB <- paste( "(", paste( paste(posnamesav, posnamesbv, sep=cb) , paste(posnamesav, negnamesbv, sep=cb), sep=sep), ")", sep="")
               COMPNAMESAB <- PNAMESAB
             })
      
    }else{
      
      if(all(negnamesb=="")){
        
        switch(method,
               "at"={
                 COMPNAMESA <- paste("(",paste(posnamesav, negnamesav, sep=sep), ")", sep="")
                 COMPNAMESAB <- paste(COMPNAMESA, posnamesbv, sep=cb)
               },
               "cellmeans"={
                 PNAMESAB <- paste( "(", paste( paste(posnamesav, posnamesbv, sep=cb) , paste(negnamesav, posnamesbv, sep=cb), sep=sep), ")", sep="")
                 COMPNAMESAB <- PNAMESAB
               })
        
      }else{
        warning("Names of contrast may be inappropriate.")
        switch(method,
               "at"={
                 COMPNAMESA <- paste("(",paste(posnamesav, negnamesav, sep=sep), ")", sep="")
                 COMPNAMESAB <- paste(paste("(",paste(COMPNAMESA, posnamesb, sep=cb), ")", sep=""), paste("(",paste(COMPNAMESA, negnamesb, sep=cb), ")", sep=""), sep=sep)
               },
               "cellmeans"={
                 PNAMESAB <- paste( "(", paste( paste(posnamesav, posnamesbv, sep=cb) , paste(negnamesav, posnamesbv, sep=cb), sep=sep), ")", sep="")
                 NNAMESAB <- paste( "(", paste( paste(posnamesav, negnamesbv, sep=cb) , paste(negnamesav, negnamesbv, sep=cb), sep=sep), ")", sep="")
                 COMPNAMESAB <- paste(PNAMESAB, NNAMESAB, sep=sep)
               })
        
      }}
    
  } 
  
  rownames(CMAB)<-COMPNAMESAB
  colnames(CMAB)<-levels(FAB)
  
  return(list(fab=FAB, cmab=CMAB))
  
}


###########################################
