## BLINK
## Installation
install.packages("devtools")

devtools::install_github("YaoZhou89/BLINK")

## Demo Data:
Demo Data could find at:  http://zzlab.net/blink/Demo_Data.zip

Note: BLINK R version only support numeric data type
## Usage
source("http://zzlab.net/GAPIT/gapit_functions.txt")

source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

library(BLINK)

### #read genotype, genotype information, phenotypes
myGM=read.table("myData.map",head=T)

myGD=read.big.matrix("myData.dat",head=F,sep="\t",type="char")

myY = read.table("myData.txt",head = T)

### # Association analysis
myBlink=Blink(Y=myY,GD=myGD,GM=myGM,maxLoop=10,time.cal=T,BIC.method="naive")
