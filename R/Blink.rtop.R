`Blink.rtop`<-function(r,df){
#Objects: transform the p value to r value
#Input: r: r value
#		df: degree of freedom, df = n-ncov-1
#		
#Output: p value
#Author: Yao Zhou
#Last Update: 8/15/2015
	tvalue=sqrt(df-1)*r/sqrt(1-r^2)
	pvalue <- 2 * pt(abs(tvalue), df-1,lower.tail = FALSE)
	return(pvalue)
}

