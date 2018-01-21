`Blink.cor`<-function(Y,GD,w=NULL,orientation="row",ms=ms,n=ny,m=nm){
  #Objects: calculate R value with covariates
  #Input: pheontype(nx1), ms is marker size for slicing the genotype, genotype(orientation="row", mxn or orientation="col", nxm,) and covariates(nxp)
  #		n is individual number, m is marker number, p is covariate number
  #Output: abs(r)
  #Author: Yao Zhou
  #Last updated: Jun 28, 2016
  if(!is.matrix(Y)) Y=as.matrix(Y)
  # Orthogonolize phenotype w.r.t. covariates
  {
    if(!is.null(w)){
      w = cbind(1,w)
    }else{
      w = matrix(1,n,1)
    }
    if(!is.matrix(w)) w = as.matrix(w)
    qw = qr(w)
    if( min(abs(diag(qr.R(qw)))) < .Machine$double.eps * m ) {
      stop("Colinear or zero covariates detected");
    }
    w = qr.Q(qw)
    tw=t(w)
    rm(qw)
  }	
  
  # Orthogonolize phenotype w.r.t. covariates
  
  {
    Y = Y - w%*%crossprod(w,Y)
    colsq = colSums(Y^2)
    div = sqrt(colsq)
    Y = Y/div
    rm(colsq,div)
  }
  time.start = proc.time()
  #Orthogonolize genotype w.r.t. covariates
  {
    if(orientation == "row"){
      rabs = matrix(NA,nrow = nrow(GD),ncol = 1)
      ntest = nrow(GD)
      ns = ceiling(ntest/ms)
      for(i in 1:ns){
        bottom=(ms*(i-1)+1)
        if(i<ns){
          up=ms*i
        }else{
          up=ntest
        }
        GDs = GD[bottom:up,]
        GDs = GDs - tcrossprod(GDs,tw)%*%tw
        rowsq= rowSums(GDs^2)
        div = sqrt(rowsq)
        GDs=GDs/div
        rabs[bottom:up,1]=abs(GDs%*%Y)
      }
    }else{
      rabs = matrix(NA,nrow = ncol(GD),ncol = 1)
      ntest = ncol(GD)
      ns = ceiling(ntest/ms)
      for(i in 1:ns){
        bottom=(ms*(i-1)+1)
        if(i<ns){
          up=ms*i
        }else{
          up=ntest
        }
        GDs = GD[,bottom:up]
        GDs = GDs- crossprod(tw,tw%*%GDs)
        colsq= colSums(GDs^2)
        div = sqrt(colsq)
        GDs=t(GDs)/div
        rabs[bottom:up,1] = abs(GDs%*%Y)
      }
    }
    rm(GDs,div)     
  }    
  time.end= proc.time()
  time.all = time.end - time.start
  print(time.all)
  return(rabs)
}
