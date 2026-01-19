SCrankcorr <- function(gamresult, computevar, ds.resolution, dsdata=FALSE){

  #### connectional rank
  ########################
  #connectional rank
  Matrix.ds<-matrix(NA, nrow=ds.resolution, ncol=ds.resolution)
  indexup.ds <- upper.tri(Matrix.ds)
  indexsave.ds <- !indexup.ds
  Matrix.ds.index<-Matrix.ds
  SClength.ds=ds.resolution*(ds.resolution+1)/2 # number of elements in the lower triangle of mat
  Matrix.ds.index[indexsave.ds]<-c(1:SClength.ds) # initial order of index in lower triangle
  Matrix.ds.SCrank<-Matrix.ds
  for (x in 1:ds.resolution){
    for (y in 1:ds.resolution){
      Matrix.ds.SCrank[x,y]<-x^2+y^2
    }
  }
  Matrix.ds.SCrank[indexup.ds]<-NA
  Matrix.ds.SCrank[indexsave.ds]<-rank(Matrix.ds.SCrank[indexsave.ds], ties.method = "average")
  gamresult.ds<-data.frame(SCrank=Matrix.ds.SCrank[indexsave.ds], computevar=NA)
  gamresult.ds$computevar <-gamresult[, computevar]
  
  correstimate <- stats::cor(gamresult.ds$SCrank, gamresult.ds$computevar, method = "spearman", use = "pairwise.complete.obs")
  p.spearman <- stats::cor.test(gamresult.ds$SCrank, gamresult.ds$computevar, method = "spearman")$p.value
  SCrankdf <- data.frame(ds.resolution=ds.resolution, Interest.var=computevar,
                         r.spearman=correstimate, p.spearman=p.spearman)
  if (dsdata == TRUE){
    names(gamresult.ds) <- c("SCrank", computevar)
    return(gamresult.ds)
  }else{
    return(SCrankdf)
  }
  
}


