# This file is part of DNAmeBloodCellAge.
# Copyright (C) 2026 Tom Bartlett and University College London
#
# DNAmeBloodCellAge is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

rm(list=ls())

set.seed(1)
libDir<-.libPaths()

CRANpackages<-c("parallel")
for(iP in 1:length(CRANpackages)){
	tempPackage<-CRANpackages[iP]
	pckg0 = try(require(tempPackage,character.only=T))
	if(!pckg0){
		pckg = try(require(tempPackage,lib.loc=libDir,character.only=T))
		if(!pckg){
			print(paste0("Installing '",tempPackage,"' from CRAN"))
			install.packages(tempPackage,lib=libDir,destdir=libDir,repos="http://cran.r-project.org")
			require(tempPackage,lib.loc=libDir,character.only=T)	
  		}
  	}
}

numCores<-detectCores()

load("BdataBloodRefCombined.rd");#Download from https://doi.org/10.5281/zenodo.20816375

Age<-as.numeric(phenData[,"Age"])/100

load("DNAmeEstPropV1.rd");#Output from DNAmeEstPropV1.R
B<-estProp[,"B"]
NK<-estProp[,"NK"]
CD4T<-estProp[,"CD4T"]
CD8T<-estProp[,"CD8T"]
Mono<-estProp[,"Mono"]
Gran<-estProp[,"Gran"]	

lmFitZ<-do.call(rbind,mclapply(1:nrow(dataB),function(tempI){	
	beta<-dataB[tempI,]			
	tempFit<-summary(lm(as.formula("beta~0+B+NK+CD4T+CD8T+Mono+Gran+Age:B+Age:NK+Age:CD4T+Age:CD8T+Age:Mono+Age:Gran")))[["coefficients"]][c("B:Age","NK:Age","CD4T:Age","CD8T:Age","Mono:Age","Gran:Age"),c("t value","Pr(>|t|)")]		
	tempZ<-as.vector((-qnorm(tempFit[,"Pr(>|t|)"]/2))*(2*(tempFit[,"t value"]>0)-1))	
	return(tempZ)
},mc.cores=numCores))
colnames(lmFitZ)<-c("B","NK","CD4T","CD8T","Mono","Gran")
rownames(lmFitZ)<-rownames(dataB)

lmFitZ[lmFitZ==Inf]<-1.1*max(lmFitZ[lmFitZ!=Inf])
lmFitZ[lmFitZ==-Inf]<-1.1*min(lmFitZ[lmFitZ!=-Inf])

save(list=c("lmFitZ","estProp","phenData"),file="DNAmeLMageV1.rd")
