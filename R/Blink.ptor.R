`Blink.ptor`<-function(p,df){
#Objects: transform the p value to r value
#Input: p: p value
#		df: degree of freedom, df = n-ncov-1
#		
#Output: r value
#Author: Yao Zhou
#Last Update: 8/15/2015
	t=qt(0.5*p,df-1,lower.tail = FALSE)
	r=sqrt(t^2/(df-1+t^2))
	return(r)
}

