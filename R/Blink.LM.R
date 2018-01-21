`Blink.LM` <-function(y,w=NULL,GDP,orientation="col"){
    N=length(y) #Total number of taxa, including missing ones
    direction=2
    if(orientation=="row"){
		GDP=t(GDP)
    }
    ntest=ncol(GDP)
    if(orientation=="row"){
        B <- matrix(NA,nrow=nrow(GDP),ncol=ncol(w)+1)
    }else{
        B <- matrix(NA,nrow=ncol(GDP),ncol=ncol(w)+1)
    }
    print(dim(B))
    if(!is.null(w)){
        nf=length(w)/N
        w=matrix(as.numeric(as.matrix(w)),N,nf  )
        w=cbind(rep(1,N),w)#add overall mean indicator
        q0=ncol(w) #Number of fixed effect excluding gnetic effects
    }else{
        w=rep(1,N)
        nf=0
        q0=1
    }
  
    y=matrix(as.numeric(as.matrix(y)),N,1  )
    
    n=N
    k=1 #number of genetic effect: 1 and 2 for A and AD respectively
    
    q1=(q0+1) # vecter index for the posistion of genetic effect (a)
    q2=(q0+1):(q0+2) # vecter index for the posistion of genetic effect (a and d)
    df=n-q0-k #residual df (this should be varied based on validating d)
    
    iXX=matrix(0,q0+k,q0+k) #Reserve the maximum size of inverse of LHS
    
    ww=crossprod(w,w)
    wy=crossprod(w,y)
    yy=crossprod(y,y)
    wwi=solve(ww)
    
    #Statistics on the reduced model without marker
    rhs=wy
    gc()
    y=as.matrix(y)
	gy=crossprod(GDP,y)
	gw=crossprod(w,GDP)
	bw=crossprod(gw,wwi)
	lbw=ncol(bw)
	P <- matrix(NA,nrow=ncol(GDP),ncol=4*(nf+1))
	for(i in 1:ntest){ 
	 	x=GDP[,i]
		xy=gy[i,1]
		xw=gw[,i]
        xx=crossprod(x,x)
     #   B21 <- crossprod(xw, wwi)
        B21=matrix(bw[i,],1,lbw)
    	t2=B21%*%xw #I have problem of using crossprod and tcrossprod here
		B22 <- xx - t2
		invB22=1/B22
		NeginvB22B21 <- crossprod(-invB22,B21)
		iXX11 <- wwi + as.numeric(invB22)*crossprod(B21,B21)
       	iXX[1:q0,1:q0]=iXX11
    	iXX[q1,q1]=invB22
        iXX[q1,1:q0]=NeginvB22B21
        iXX[1:q0,q1]=NeginvB22B21
        rhs=c(wy,xy)
        beta <- crossprod(iXX,rhs)   #both a and d go in
        df=n-q0-1
        ve=(yy-crossprod(beta,rhs))/df
        se=sqrt(diag(iXX)*ve)
        tvalue=beta/se
        pvalue <- 2 * pt(abs(tvalue), df,lower.tail = FALSE)
        if(!is.na(abs(B22[1,1]))){
            if(abs(B22[1,1])<10e-8)pvalue[]=NA}
        P[i,]=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
        print(length(beta))
    #   	B[i,]=beta[length(beta)]
    }
    return(list(P=P,PF=NULL,beta=beta))
} #end of function


