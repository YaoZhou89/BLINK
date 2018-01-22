# How to run multiple traits?
### BLINKR could only run one trait per time, if you want to run all trails in the myY, you can run like this:
    for (i in 2:ncol(myY)){

      myBlink=Blink(Y=myY[,c(1,i)],GD=myGD,GM=myGM,maxLoop=10,time.cal=T,BIC.method="naive")
  
    }

# How to add coviates?
### Note: the covriates file should be n by q, n is the individuals number and q is the number of covriates, and the order should be the same with myY and myGD. Like adding the top 3 PCs as coviates:
### 1. calculate PCs using prcomp:
    cov = prcomp(as.matrix(myGD))
### 2. then extract the top 3 PCs:
    myCV = cov$x[,1:3]
### 3. run BLINK with PCs:
    myBlink=Blink(Y=myY,GD=myGD,GM=myGM,CV = myCV, maxLoop = 10, time.cal = T)
    
## Parameters may be useful
    time.cal: calculate the time spend in each step;

    maxLoop: the loop for iteration;

    maf.threshold: filter SNPs with MAF < maf.threshold
