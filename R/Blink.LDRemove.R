`Blink.LDRemoveBlock`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
#`Blink.LDRemove`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
#Objects: Calculate LD and remove the correlated SNPs
#Authors: Yao Zhou
#Last Update:  03/03/16
	if (model=="D"){
		GDneo=1-abs(GDneo-1)
	}

	GDneo=as.matrix(GDneo)
	if(min(ncol(GDneo),nrow(GDneo))<201) bound=FALSE
	if(orientation=="col"){
		n=nrow(GDneo)
		if(bound){
			GDneo=GDneo[sample(n,200,replace=F),]
		}
	}else{
		n=ncol(GDneo)
		if(bound){
			GDneo=GDneo[,sample(n,200,replace=F)]
		}
		GDneo=t(GDneo)
	}
	# cat("ncol(GDneo) is",ncol(GDneo),"\n")
	corr=cor(GDneo)
	corr[is.na(corr)]=1
	corr[abs(corr)<=LD]=0
	corr[abs(corr)>LD]=1
	Psort=as.numeric(matrix(1,1,ncol(corr)))
	# print(ncol(corr))
	for(i in 2:ncol(corr)){
		p.a=Psort[1:(i-1)]
		p.b=as.numeric(corr[1:(i-1),i])
		index=(p.a==p.b)
		index[(p.a==0)&(p.b==0)]=FALSE
		if(sum(index)!=0) Psort[i]=0
	}
	seqQTN=Porder[Psort==1]
	return(seqQTN)
}

`Blink.LDRemove`<-function(GDneo=NULL,LD=0.7,Porder=NULL,bound=FALSE,model="A",orientation="row",block=1000,LD.num =50){
#Objects: LD remove, especially length(Porder)>10000
#Authors: Yao Zhou
#Last update: 08/15/2016
  #GDneo = as.matrix(GDneo)
  GDneo = as.matrix(as.data.frame(GDneo))
  SNP.index = apply(GDneo,1,sd)!=0
  GDneo = GDneo[SNP.index,]
  Porder = Porder[SNP.index]
  l = block
	seqQTN=NULL
	lp=length(Porder)
	k=ceiling(lp/l)
	GDneo=as.matrix(GDneo)
	if(min(ncol(GDneo),nrow(GDneo))<201) bound=FALSE
	if(orientation=="col"){
		n=nrow(GDneo)
		if(bound){
			GDneo=GDneo[sample(n,200,replace=F),]
		}
	}else{
		n=ncol(GDneo)
		if(bound){
			GDneo=GDneo[,sample(n,200,replace=F)]
		}
		GDneo=t(GDneo)
	}
	for(i in 1:k){
		bottom=(i-1)*l+1
		up=l*i
		if(up>lp) up = lp
		Porderb=Porder[bottom:up]

		index = seq(bottom:up)
		GDneob = GDneo[,index]
		# cat("i is ",i,"\n")
		# print(length(index))
		seqQTNs = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb,orientation="col",model=model)
		# print(seqQTN)
		seqQTN = append(seqQTN,seqQTNs)
		if(k >1){
		  index1 = which(Porder %in% seqQTN)
		  Porderb = Porder[index1]
		  GDneob = GDneo[,index1]
		  if(length(index1)>1){
		    seqQTN = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb,orientation="col",model=model)
		  }else{
		    seqQTN = Porderb
		  }

		}
		if(LD.num < length(seqQTN)) break
	}
	rm(GDneob,Porderb)
	return(seqQTN)
}
#
# `Blink.LDRemove`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
# #`Blink.LDRemovebackup`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
# #Objects: LD remove, especially length(Porder)>10000
# #Authors: Yao Zhou
# #Last update: 08/15/2016
# 	seqQTN=NULL
# 	is.done=FALSE
# 	l=1000
# 	lp=length(Porder)
# 	tt=1
# 	n = ncol(GDneo)
# 	n.limit = n/log(n)
# 	while(!is.done){
# 		tt = tt+1
# 		Pordern=Blink.LDRemoveDivided(GDneo=GDneo,LD=LD,Porder=Porder,orientation=orientation,model=model,l=l)
# 		index=Porder %in% Pordern
# 		if(orientation=="col"){
# 			GDneo = GDneo[,index]
# 		}else{
# 			GDneo = GDneo[index,]
# 		}
# 		ls=length(Pordern)
# 		if(ls==lp) lp=l*tt
# 		if(ls<=lp){
# 			is.done=TRUE
# 		}
# 		Porder = Pordern
# 	}
# 	if(length(Porder) > 1){
# 		seqQTN=Blink.LDRemoveBlock(GDneo=GDneo,LD=LD,Porder=Porder,orientation=orientation,model=model)
# 	}else{
# 		seqQTN = Porder
# 	}
# 	return(seqQTN)
# }

