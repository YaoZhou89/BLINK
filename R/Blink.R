`Blink` <- function(Y=NULL,QTN.position=NULL,GD=NULL,GM=NULL,CV=NULL,DPP=100000000,kinship.algorithm="FARM-CPU",file.output=TRUE,cutOff=0.01,method.GLM="FarmCPU.LM",Prior=NULL,ncpus=1,maxLoop=10,LD=0.7,threshold.output=.0001,alpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),WS=c(1e0,1e3,1e4,1e5,1e6,1e7),GP=NULL
,maxOut=10,converge=1,iteration.output=FALSE,acceleration=0,threshold=NA,model="A",MAF.calculate=FALSE,plot.style="FarmCPU",p.threshold=NA,maf.threshold=0,bound=FALSE,method.sub="reward",method.sub.final="reward",stepwise=FALSE,BIC.method="even",LD.wise=FALSE,time.cal=FALSE,Prediction = F){
	print("----------------------Welcome to Blink----------------------")
	time.start=proc.time()
	nm=nrow(GM)
	ny=nrow(Y)
	if(is.na(threshold)) threshold=floor(ny/log(ny))
	ngd=nrow(GD)
	orientation="col"
	theSNP=2
	ns=nrow(GD)
	seqQTN=NULL
	if(nm==ngd){
		orientation="row"
		theSNP=1
		ns=ncol(GD)
	}
	if(maf.threshold > 0) {
		MAF.calculate = TRUE
	}
	
	if(MAF.calculate==FALSE){
        MAF=NA
    }else{
        MAF=apply(GD,theSNP,mean)
        MAF=matrix(MAF,nrow=1)
        MAF=apply(MAF,2,function(x) min(1-x/2,x/2))
    }
    MAF.index = 1:nm
    if(maf.threshold > 0) {
		MAF.index = MAF > maf.threshold	
		MAF = MAF[MAF.index]
	}
	ac=NULL
	if(acceleration!=0) ac=rep(1.0,nm)
	index=which(GM[,3]==0 )
	if(length(index)>0){
		GM[index,3]=1
	 }

	P=GP
	gc()
	if(ncol(GD)>nm & orientation=="col"){
		if(is.big.matrix(GD)){
			GD=deepcopy(GD,rows=1:nrow(GD),cols=2:ncol(GD))
		}else{
			GD=as.matrix(GD[,-1])
		}
	}
	# GD=as.matrix(GD)
	gc()
	shift=0
	for(trait in 2:ncol(Y)){
		name.of.trait=colnames(Y)[trait]
		index=MAF.index
		seqTaxa=which(!is.na(Y[,trait]))
		Y1=Y[seqTaxa,]
		if(!is.null(CV)){
			if(ncol(CV)==1){
				CV1=CV[seqTaxa]
			}else{
				CV1=CV[seqTaxa,]
			}
		}else{
			CV1=NULL
		}
		
		if(orientation=="col"){
            if(is.big.matrix(GD)){
                GD1=deepcopy(GD,rows=seqTaxa,cols=index)
            }else{
                GD1=GD[seqTaxa,index]
            }
        }else{
            if(is.big.matrix(GD)){
                GD1=deepcopy(GD,rows=index,cols=seqTaxa)
            }else{
                GD1=GD[index,seqTaxa]
                GD1=as.matrix(GD1)
            }
        }
		LD.time=rep(0,maxLoop)
		BIC.time=rep(0,maxLoop)
		GLM.time=rep(0,maxLoop)
		theLoop=0
		theConverge=0
		seqQTN.save=c(0)
		isDone=FALSE
		name.of.trait2=name.of.trait
		while(!isDone) {
			theLoop=theLoop+1
			print(paste("----------------------Iteration:",theLoop,"----------------------",sep=" "))
			if(iteration.output) name.of.trait2=paste("Iteration_",theLoop,".",name.of.trait,sep="")
			myPrior=FarmCPU.Prior(GM=GM,P=P,Prior=Prior,kinship.algorithm=kinship.algorithm)	
   			if(!is.null(myPrior)){
				if(theLoop!=1){
  	  				seqQTN.p=myPrior
  	      			if(theLoop==2){
        				index.p=seqQTN.p<(0.01/nm)
            			if(!is.na(p.threshold)){
               				index.p=seqQTN.p<p.threshold
            			}
						index.p[is.na(index.p)]=FALSE
            			seqQTN.selected=as.numeric(which(index.p))
        			}else{
          	  			index.p=seqQTN.p<(1/nm)
          	  			if(!is.na(p.threshold)){
               				index.p=seqQTN.p<p.threshold
            			}
						index.p[is.na(index.p)]=FALSE
            			seqQTN.selected=as.numeric(which(index.p))		
        			}

					Porder=order(myPrior[seqQTN.selected],na.last=T,decreasing=FALSE)
					t1=proc.time()
					if(length(Porder)>1){
						if(LD.wise & (length(Porder)>threshold)){
							max_Porder=max(Porder)
							if(max_Porder > 10000) max_Porder=10000
							step_bin= 10
							Porder_new=rep(Porder[1],threshold)
							Po=1
							for( po in 2:max_Porder){
								if(min(abs(seqQTN.selected[Porder[po]]-seqQTN.selected[Porder_new]))>step_bin){
									Po=Po+1
									Porder_new[Po]=Porder[po]
									if (Po >=threshold) break
								}
							}
							Porder=Porder_new
            			}
            			
						if(is.big.matrix(GD1)){
							if(orientation=="col"){
								GDnew=deepcopy(GD1,cols=seqQTN.selected)
								GDneo=deepcopy(GDnew,cols=Porder)
							}else{
								GDnew=deepcopy(GD1,rows=seqQTN.selected)
								GDneo=deepcopy(GDnew,rows=Porder)						
							}
						}else{
							if(orientation=="col"){
								GDnew=GD1[,seqQTN.selected]					
								GDneo=GDnew[,Porder]											
							}else{
								GDnew=GD1[seqQTN.selected,]					
								GDneo=GDnew[Porder,]						
							}		
						}
						
						print("LD remove is working....")
						print("Number SNPs for LD remove:")
						print(length(Porder))
						Psort=Blink.LDRemove(Porder=Porder,GDneo=GDneo,bound=bound,LD=LD,model=model,orientation=orientation)
						seqQTN.can=seqQTN.selected[Psort]
						t2=proc.time()
						print("Model selection based on BIC is working....")
						print("Number of SNPs for BIC selection:")
						print(length(seqQTN.can))
						myBIC = Blink.BICselection(Psort=seqQTN.can,GD=GD1,Y=Y1,orientation=orientation,BIC.method=BIC.method)
						seqQTN = myBIC$seqQTN
						#if(theLoop==6) print(seqQTN)
						t3=proc.time()
						LD.time[theLoop]=as.numeric(t2)[3]-as.numeric(t1)[3]
						BIC.time[theLoop]=as.numeric(t3)[3]-as.numeric(t2)[3]
					}else if(length(Porder)==1){
						print("LD remove is working....")
						print("Model selection based on BIC is working....")
						seqQTN=seqQTN.selected
					}else{
						seqQTN=NULL
					}
				}
			}else{
				seqQTN=NULL
			}
		
			if(theLoop==2){
       			if(!is.na(p.threshold)){
        			if(min(myPrior,na.rm=TRUE)>p.threshold){
         	    		seqQTN=NULL
          		   		print("Top snps have little effect, set seqQTN to NULL!")
          			}
        		}else{
        			if(min(myPrior,na.rm=TRUE)>0.01/nm){
						seqQTN=NULL
       		     		print("Top snps have little effect, set seqQTN to NULL!")
       		 		}
          	 	}
    		}	
			if(theLoop==2&&is.null(seqQTN)|length(seqQTN)==0&&theLoop==2){
	   			print(paste("seqQTN is:",seqQTN,",stop here",sep=""))
    			if(!isDone | iteration.output){
					gc()
        			p.GLM=GWAS[,4]
            		p.GLM.log=-log10(quantile(p.GLM,na.rm=TRUE,0.05))
            		bonf.log=1.3
            		bonf.compare=p.GLM.log/bonf.log
            		p.FARMCPU.log=-log10(p.GLM)/bonf.compare
        			GWAS[,4]=10^(-p.FARMCPU.log)
            		GWAS[,4][which(GWAS[,4]>1)]=1
            		colnames(GWAS)=c(colnames(GM),"P.value","maf","nobs","Rsquare.of.Model.without.SNP","Rsquare.of.Model.with.SNP","FDR_Adjusted_P-values")
            		Vp=var(Y1[,2],na.rm=TRUE)
            		if(file.output) GAPIT.Report(name.of.trait=name.of.trait2,GWAS=GWAS,pred=NULL,tvalue=NULL,stderr=stderr,Vp=Vp,DPP=DPP,cutOff=cutOff,threshold.output=threshold.output,MAF=MAF,seqQTN=QTN.position,MAF.calculate=MAF.calculate,plot.style=plot.style)
                    myPower=GAPIT.Power(WS=WS, alpha=alpha, maxOut=maxOut,seqQTN=QTN.position,GM=GM,GWAS=GWAS,MaxBP=1e10)
        		} 
            	break
			}
			if(theLoop>1){
				if(seqQTN.save!=0 & seqQTN.save!=-1 & !is.null(seqQTN)) seqQTN=union(seqQTN,seqQTN.save)
			}
			if(theLoop>2 ){
				if( length(Porder)>1){
					BIC=Blink.BICselection(Psort=seqQTN,GD=GD1,Y=Y1,orientation=orientation,BIC.method=BIC.method)
					seqQTN = BIC$seqQTN
				}
			}
	
			print("seqQTN:")   
			print(seqQTN)
			theConverge=length(intersect(seqQTN,seqQTN.save))/length(union(seqQTN,seqQTN.save))
			isDone=((theLoop>=maxLoop)|(theConverge>=converge))
			if(!is.null(seqQTN)) seqQTN.save=seqQTN
			gc()
			if(!is.null(seqQTN)){
				if(orientation=="col"){
					theCV=cbind(CV1,GD1[,seqQTN])
				}else{
					if(length(seqQTN)>1){
						theCV1=t(GD1[seqQTN,])
						theCV=cbind(CV1,theCV1)
					}else{
						theCV=cbind(CV1,GD1[seqQTN,])
					}
				}
			}else{
				theCV=CV1
			}
			t4=proc.time()
			if(!is.null(theCV)) theCV=as.matrix(theCV)
			#if(theLoop==4) write.table(theCV,"CV.txt",col.names=F,row.names=F,quote=F,sep="\t")
			myGLM=FarmCPU.LM(y=Y1[,trait],GDP=GD1,w=theCV,orientation=orientation)
   			#print(dim(myGLM$P))
		    #myGLM=Blink.LM(y=Y1[,trait],GDP=GD1,w=theCV,orientation=orientation)
   			#print(dim(myGLM$P))

   			if(!isDone){
    			myGLM=Blink.SUB(GM=GM,GLM=myGLM,QTN=GM[seqQTN,],method=method.sub,model=model)
    		}else{
    			#save(myGLM,file="myGLM_last.Rdata")
        		myGLM=Blink.SUB(GM=GM,GLM=myGLM,QTN=GM[seqQTN,],method=method.sub,model=model)
    		}
			t5=proc.time()
			GLM.time[theLoop]=as.numeric(t5)[3]-as.numeric(t4)[3]
			P=myGLM$P[,ncol(myGLM$P)-shift]
			index=which(ac>1)
			P[P==0] <- min(P[P!=0],na.rm=TRUE)*0.01
			P[is.na(P)] =1
			gc()
			nf=ncol(myGLM$P)/4
			tvalue=myGLM$P[,nf*2-shift]
			stderr=myGLM$P[,3*nf-shift]
			GWAS=cbind(GM[MAF.index,],P,MAF,NA,NA,NA,NA)
			colnames(GWAS)=c(colnames(GM),"P.value","maf","nobs","Rsquare.of.Model.without.SNP","Rsquare.of.Model.with.SNP","FDR_Adjusted_P-values")
			Vp=var(Y1[,2],na.rm=TRUE)		
			if(file.output){
				if(theLoop==1&&is.null(CV)){
					GAPIT.Report(name.of.trait=name.of.trait2,GWAS=GWAS,pred=NULL,ypred=NULL,tvalue=NULL,stderr=stderr,Vp=Vp,DPP=DPP,cutOff=cutOff,threshold.output=threshold.output,MAF=MAF,seqQTN=QTN.position,MAF.calculate=MAF.calculate,plot.style=plot.style)
				}else{
					GAPIT.Report(name.of.trait=name.of.trait2,GWAS=GWAS,pred=NULL,ypred=NULL,tvalue=NULL,stderr=stderr,Vp=Vp,DPP=DPP,cutOff=cutOff,threshold.output=threshold.output,MAF=MAF,seqQTN=QTN.position,MAF.calculate=MAF.calculate,plot.style=plot.style)
				}
			}# end of file.out		
		}	#end of theLoop
		PEV=NULL
		if(Prediction){
			YP = cbind(Y[,1],Y[,trait])
			p.rank = order(GWAS[,4],na.last = T,decreasing=F)
			if(sum(GWAS[p.rank,4]<p.threshold)>0){
				seqQTN = p.rank[GWAS[p.rank,4] < p.threshold]
			}else{
				seqQTN = p.rank[1]
			}
			if(!is.big.matrix(GD)){
				if(orientation=="col"){
					GDpred = GD[,seqQTN]
				}else{
					if(length(seqQTN)>1){
						GDpred = t(GD[seqQTN,])
					}else{
						GDpred = as.matrix(GD[seqQTN,])
					}
				}
			}else{
				if(orientation=="col"){
					GDpred = deepcopy(GD,cols=seqQTN)
				}else{
					GDpred = deepcopy(GD,rows=seqQTN)
				}
			}
			PEV = Blink.Pred(Y = YP,GD = GDpred, CV = CV,orientation = orientation)
		}
		if(time.cal){
			print("LD.time(sec):")
			print(LD.time[1:theLoop])
			print("BIC.time(sec):")
			print(BIC.time[1:theLoop])
			print("GLM.time(sec):")
			print(GLM.time[1:theLoop])
		}
		time.end=proc.time()
		time.all=as.numeric(time.end)[3]-as.numeric(time.start)[3]
		print(paste("-------------Blink finished successfully in",round(time.all,2),"seconds!-----------------"))
	#	print(proc.time())
		write.table(GWAS,paste(name.of.trait2,"_GWAS.txt",sep=""),sep="\t",col.names=T,row.names=F)
	}#end of phenotype
	return(list(GWAS=GWAS,myGLM=myGLM,PEV = PEV))
}#	end of function Blink	

