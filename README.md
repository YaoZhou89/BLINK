## BLINK
This is the R version of BLINK model. You can find the C version [here](https://github.com/Menggg/BLINK)
## Installation
    install.packages("devtools")

    devtools::install_github("YaoZhou89/BLINK")

## Demo Data
BLINK R version only support numeric data type.

## Usage
#### # source functions needed
    source("http://zzlab.net/GAPIT/gapit_functions.txt")

    source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

    library(BLINK)

#### # genotype information data
    myGM=read.table("myData.map",head=T)
#### # genotype data
    myGD=read.big.matrix("myData.dat",head=F,sep="\t",type="char") 
#### # phenotype data
    myY = read.table("myData.txt",head = T) 

#### # association analysis
    myBlink=Blink(Y=myY,GD=myGD,GM=myGM,maxLoop=10,time.cal=T) 
#### # more features   
More parameters explained [here](https://github.com/YaoZhou89/BLINK/blob/master/man/User%20Manual.md)

## The license notice
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


## Author
Dr. Yao Zhou (yao.zhou@genetics.ac.cn)
