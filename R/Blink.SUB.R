`Blink.SUB` <-
function(GM=NULL,GLM=NULL,QTN=NULL,method="mean",useapply=TRUE,model="A"){
    #Input: FarmCPU.GLM object
    #Input: QTN - s by 3  matrix for SNP name, chromosome and BP
    #Input: method - options are "penalty", "reward","mean","median",and "onsite"
    #Requirement: P has row name of SNP. s<=t. covariates of QTNs are next to SNP
    #Output: GLM with the last column of P updated by the substituded p values
    #Authors: Xiaolei Liu and Zhiwu Zhang
    # Last update: Febuary 26, 2013
    ##############################################################################
    if(is.null(GLM$P)) return(NULL)  #P is required
    if(is.null(QTN)) return(NULL)  #QTN is required
    position=match(QTN[,1], GM[,1], nomatch = 0)
    nqtn=length(position)
    if(model=="A"){
        index=(ncol(GLM$P)-nqtn):(ncol(GLM$P)-1)
        spot=ncol(GLM$P)
    }else{
        index=(ncol(GLM$P)-nqtn-1):(ncol(GLM$P)-2)
        spot=ncol(GLM$P)-1
    }
    if(ncol(GLM$P)!=1){
        if(length(index)>1){
            if(method=="penalty") P.QTN=apply(GLM$P[,index],2,max,na.rm=TRUE)
            if(method=="reward") P.QTN=apply(GLM$P[,index],2,min,na.rm=TRUE)
            if(method=="mean") P.QTN=apply(GLM$P[,index],2,mean,na.rm=TRUE)
            if(method=="median") P.QTN=apply(GLM$P[,index],2,median,na.rm=TRUE)
            if(method=="onsite") P.QTN=GLM$P0[(length(GLM$P0)-nqtn+1):length(GLM$P0)]
        }else{
            if(method=="penalty") P.QTN=max(GLM$P[,index],na.rm=TRUE)
            if(method=="reward") P.QTN=min(GLM$P[,index],na.rm=TRUE)
            if(method=="mean") P.QTN=mean(GLM$P[,index],na.rm=TRUE)
            if(method=="median") P.QTN=median(GLM$P[,index],median,na.rm=TRUE)
            if(method=="onsite") P.QTN=GLM$P0[(length(GLM$P0)-nqtn+1):length(GLM$P0)]
        }
        GLM$P[position,spot]=P.QTN
    }
    return(GLM)
}#The function FarmCPU.SUB ends here

