# if (!require(softImpute)) install.packages('softImpute')
# if (!require(missForest)) install.packages('missForest')
# if (!require(missMDA)) install.packages('missMDA')
# if (!require(bnstruct)) install.packages('bnstruct')
#

library(softImpute)
library(missForest)
library(missMDA)
library(bnstruct)
#

### master programme
missHyb <- function(X, p1 = 'MF', p2 = 'MDA', mat.init = NULL,n.pca = NULL,k.nn=15,cat.var = c()) {

  ### phase 1 methods
  if(p1 == 'MF') {
    ximp1 = missForest(X)$ximp
  } else if(p1 == 'MDA') {
    if(is.null(n.pca)){
      n.pca <- estim_ncpPCA(X, ncp.max=30)$ncp
    }
      ximp1 = imputePCA(X, ncp=n.pca)$completeObs
  } else if (p1 == 'KNN'){
      ximp1 = knn.impute(X, k= k.nn, cat.var = cat.var )
  } else if(p1 == 'SI'){
    soft.rank<-c(30)
    #lam0<-lambda0(cbind(X))
    #lam0
    soft.lambda<-c(1)
    soft.g11<-softImpute(cbind(X),rank.max=soft.rank[1],lambda=soft.lambda[1],trace=FALSE,type="svd")
    #soft.pred11<-softImpute::complete(cbind(ZM110),soft.g11)*rep.row(sd11,dim(ZM110)[1])+rep.row(mean11,dim(ZM110)[1])
    ximp1 = softImpute::complete(cbind(X),soft.g11)
  }


  #### phase 2 methods:
  if(p2 == 'MF') {
    ximp = impute_MDA_MF(dat, mat.init = ximp1)
  }
  else if(p2 == 'MDA') {
    # if(is.null(n.pca)){
    #   n.pca <- estim_ncpPCA(X, ncp.max=30)
    # }
    ximp = impute_MF_MDA(dat, mat.init =  ximp1, n.pca = n.pca)
  }
  else if (p2 == 'KNN'){
    ximp = impute_MDA_KNN(dat, mat.init = ximp1, k = k.nn)
  }
  else if(p2 == 'SI'){
    ximp = impute_MF_SI(dat, mat.init = ximp1)
  }

  return(ximp)
}



### first missForest and then followed by missMDA
impute_MF_MDA <- function(ZM110, mat.init = NULL, n.pca=NULL){
  #p <- dim(ZM110)[2]

  #Yimp <- mat.init

  if(is.null(n.pca)){
    n.pca <- estim_ncpPCA(ZM110,ncp.max=30)$ncp
  }

  idx_miss = which(apply(is.na(ZM110),2,sum) > 0)

  for(j in idx_miss){
    mat.tmp = mat.init
    mat.tmp[,j] = ZM110[,j]
    mda.g11.j <- imputePCA(mat.tmp, ncp=n.pca)$completeObs
    mat.init[,j] = mda.g11.j[,j]
  }
  return(mat.init)

}


impute_MF_SI <- function(ZM110, mat.init = NULL){

  idx_miss = which(apply(is.na(ZM110),2,sum) > 0)

  for(j in idx_miss){
    mat.tmp = mat.init
    mat.tmp[,j] = ZM110[,j]
    #mda.g11.j <- imputePCA(mat.tmp, ncp=n.pca)$completeObs

    soft.rank<-c(20)
    lam0<-lambda0(cbind(ZM110))
    #lam0
    soft.lambda<-c(1)
    soft.g11<-softImpute(cbind(mat.tmp),rank.max=soft.rank[1],lambda=soft.lambda[1],trace=FALSE,type="svd")
    #soft.pred11<-softImpute::complete(cbind(ZM110),soft.g11)*rep.row(sd11,dim(ZM110)[1])+rep.row(mean11,dim(ZM110)[1])
    soft.g11 = softImpute::complete(cbind(mat.tmp),soft.g11)
    mat.init[,j] = soft.g11[,j]


    #mat.init[,j] = mda.g11.j[,j]
  }
  return(mat.init)

}




### first missMDA and then use missForest,
### the mat.init use the output from missMDA or other method
### takes too long, not practical
impute_MDA_MF <- function(ZM110, mat.init = NULL, n.pca=NULL){
  #p <- dim(ZM110)[2]

  #Yimp <- mat.init

  # if(is.null(n.pca)){
  #   n.pca <- estim_ncpPCA(ZM110,ncp.max=30)
  # }

  idx_miss = which(apply(is.na(ZM110),2,sum) > 0)

  for(j in idx_miss){
    t0 = Sys.time()
    mat.tmp = mat.init
    mat.tmp[,j] = ZM110[,j]
    #mda.g11.j <- imputePCA(mat.tmp, ncp=n.pca)$completeObs
    mda.g11.j <- missForest::missForest(xmis=mat.tmp, maxiter = 10)$ximp
    mat.init[,j] = mda.g11.j[,j]
    #print(paste(j,Sys.time() - t0))
  }
  return(mat.init)

}



### first missMDA and then use KNN,
### the mat.init use the output from missMDA or other method
### takes too long, not practical
impute_MDA_KNN <- function(ZM110, mat.init = NULL, n.pca=NULL,k.nn=15){
  #p <- dim(ZM110)[2]

  #Yimp <- mat.init

  # if(is.null(n.pca)){
  #   n.pca <- estim_ncpPCA(ZM110,ncp.max=30)
  # }

  idx_miss = which(apply(is.na(ZM110),2,sum) > 0)

  for(j in idx_miss){
    t0 = Sys.time()
    mat.tmp = mat.init
    mat.tmp[,j] = ZM110[,j]
    #mda.g11.j <- imputePCA(mat.tmp, ncp=n.pca)$completeObs
    #mda.g11.j <- missForest::missForest(xmis=ZM110, maxiter = 1)$ximp
    mda.g11.j <- knn.impute(mat.tmp, k= k.nn, cat.var=c() )
    mat.init[,j] = mda.g11.j[,j]
    #print(paste(j,Sys.time() - t0))
  }
  return(mat.init)

}
