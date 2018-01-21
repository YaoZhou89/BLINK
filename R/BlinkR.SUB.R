`BlinkR.SUB` <-function(CV,seqQTN,Y,r,ny,ms,m){
#Objects: subsitution of r value for covariates
#Input: Y, nx1 vector, phenotype
#		w, all covariates without 1 
#		seq = length(seqQTN), number of SNPs added as covariate
#Output：r, r value for SNPs added as covariates
#Author:Yao Zhou
#Last update: 08/15/2016
	rsnp=matrix(NA,nsnp,1)
	ncov=ncol(CV)
	nf=ncov-nsnp
	GDP=CV[,(nf+1):ncov]
	for(i in 1:nsnp){
		w=GDP[,-i]
		if(nsnp==1) w=NULL
		GD=as.matrix(GDP[,i])
		rsnp[i,1]=Blink.cor(Y=Y,w=w,GD=GD,orientation=orientation,,ms=ms,n=ny,m=nm)
	}
	rm(GDP, GD, w, ncov, nf)
    return(rsnp)
}