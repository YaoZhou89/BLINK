`Blink.BICselection` <-  function(Psort=NULL,CV=NULL,GD=NULL,Y=Y1,orientation=NULL,BIC.method="even"){
#Objects: fixed model selection using BIC
#Input:Y,GD,Psort
#		BIC.method: Naive: detect all SNPs of Psort
#					even: detect some SNPs, step=floor(sqrt(m))+1
#					fixed: detect the SNPs by fixed steps. Default is 20
#					log: detect the SNPs by log(10,N) transform
#					ln: detect the SNPs by ln transform
#Output: seqQTN: SNP position
#Author: Yao Zhou
#Last update: 01/05/2016, modified 03/31/2016
	GD = as.matrix(GD)
	n=nrow(Y)
	threshold=floor(n/log(n))
	if(threshold < length(Psort)){
		seqQTN=Psort[1:threshold]
	}else{
		seqQTN=Psort
	}
	y=Y[,2]

	s=0
	a=0
	m=length(seqQTN)
	pmatrix=matrix(1,m,m)

	if(BIC.method=="naive"){
		position=seq(1:m)

	} else if(BIC.method=="even"){
		step.length=floor(sqrt(m))+1
		step=floor(m/step.length)
		if ((m-step*step.length)>=(0.5*step.length)) {
        	step=step+1
   		}
		if (step.length>m) {
			step.length=m
      	  	step=1
   		}
		position=seq(step,m,step)
		if(position[length(position)]<m) position=append(position,m)
	} else if(BIC.method=="lg"){
		if(m==1){
			position =c(1)
		}else{
			le=seq(1:m)
			step=le/log10(le)
			for(i in 2:m){
				le[i]=le[i-1]+step[i]
				le=round(le)
				if(le[i]>m){
					position=le[1:i]
					break
				}
			}
		}

	} else if(BIC.method=="ln"){
		if(m==1){
			position =c(1)
		}else{
			le=seq(1:m)
			step=le/log(le)

			for(i in 2:m){
				le[i]=le[i-1]+step[i]
				le=round(le)
				if(le[i]>m){
					position=le[1:i]
					break
				}
			}
		}
	} else if(BIC.method=="fixed"){
		if(m>20){
			position=floor(seq(1,m,m/20))
		}else{
			position=seq(1:m)
		}
	}else{
		print("please choose one method for BIC")
		break
	}

	BICv=rep(NA,length(position))
	if(is.null(CV)){
		w=as.matrix(rep(1,n))
		ww=n
		ncov=2
	}else{
		CV=as.matrix(CV)
		w=cbind(1,CV)
		ww=crossprod(w)
		ncov=ncol(ww)+1
	}
 	wwi=ginv(ww)


	pos.pre=0
	k=0
	for(pos in position){
		if(pos>m) pos=m
		if(orientation=="col"){
			x=GD[,seqQTN[(pos.pre+1):pos]]
		}else{
			x=GD[seqQTN[(pos.pre+1):pos],]
			if(is.matrix(x)){
				x=t(x)
			}else{
				x=as.matrix(x)
			}
		}
		k=k+1
		pos.pre=pos
		x=as.matrix(x)

		if(k==1){
			ww=crossprod(w,w)
		}else{
			WW=matrix(0,(nwc+nxc),(nwc+nxc))
			WW[1:nwc,1:nwc]=ww
			WW[1:nwc,(nwc+1):(nwc+nxc)]=xw
			WW[(nwc+1):(nxc+nwc),1:nwc]=wx
			WW[(nwc+1):(nwc+nxc),(nwc+1):(nwc+nxc)]=xx
			ww=WW
		}
		nwc = ncol(w)
		nxc = ncol(x)
  	  	iXX = matrix(0,(nwc+nxc),(nwc+nxc))
		xx = crossprod(x,x)
		xw = crossprod(x,w)
		wx = crossprod(w,x)
		t1 = wwi %*% wx
		t2 = xx - xw %*% t1
		if (!is.null(t2)){
		M22 = ginv(t2)
		t3=xw %*% wwi
		M21=-M22 %*% t3
		M12=-t1 %*% M22
		M11=wwi + t1 %*% M22 %*% t3
		iXX[1:nwc,1:nwc]=M11
		iXX[(nwc+1):(nwc+nxc),(nwc+1):(nwc+nxc)]=M22
		iXX[(nwc+1):(nwc+nxc),1:nwc]=M21
		iXX[1:nwc,(nwc+1):(nwc+nxc)]=M12
		w=cbind(w,x)
		wy=crossprod(w,y)
		wwi=iXX
		beta=wwi %*% wy
		yp= w %*% beta
		ve=as.numeric(var(yp-y))
		RSS= (yp-y)^2
		n2LL=n*log(2*pi)+n*log(ve)+2*sum(RSS/(2*ve))
		# BICv[k]=n2LL+2*(nwc+nxc-1)*log(n)
		BICv[k]=n2LL+(nwc+nxc-1)*log(n)
		df=(n-pos-1)
		MSE=sum(RSS)/df
		se=sqrt(diag(iXX)*MSE)
		tvalue=beta/se
        pvalue <- 2 * pt(abs(tvalue), df,lower.tail = FALSE)
        pmatrix[1:pos,pos]=pvalue[ncov:length(pvalue)]
    }
	}
	seqQTN=Psort[1:position[which(BICv==min(BICv,na.rm=T))]]
	pvalue=as.numeric(pmatrix[1:length(seqQTN),length(seqQTN)])
	return(list(seqQTN=seqQTN,pvalue=pvalue,BIC=BICv))
}


