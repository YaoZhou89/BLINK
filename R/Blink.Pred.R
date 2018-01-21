`Blink.Pred` <- function(GD = NULL, Y = NULL,CV = NULL){
## Objects: Prediction using significant pseudo QTNs
## Input: Y, CV and GD
## Output: Predicted Phenotype
## Authors: Yao Zhou
## Last update: 2/6/2017

	if(is.big.matrix(GD)) GD = as.matrix(GD)
	if(orientation =="row"){
		GD = t(GD)
		if(nrow(GD)==1) GD = t(GD)
	}
	
	seqTaxa=which(!is.na(Y[,2]))
	Y1=Y[seqTaxa,2]
	GD1 = GD[seqTaxa,]
	
	if(is.null(CV)){
		mylm = lm(Y1 ~ GD1)
		PEV = predict(mylm,as.data.frame(GD))
	}else{
		CV1 = CV[seqTaxa,]
		mylm = lm(Y1 ~ CV1 + GD1)
		PEV = predict(mylm,as.data.frame(cbind(CV,GD)))
	}
	return(PEV)	
}

