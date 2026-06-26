# This file is part of DNAmeBloodCellAge.
# Copyright (C) 2026 University College London.
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
nBoot<-10000
nFeat<-1000
libDir<-.libPaths()

CRANpackages<-c("parallel","MCMCpack")
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

load("DNAmeLMageV1.rd");#Output from DNAmeLMageV1.R
topFeatures<-unique(as.vector(unlist(apply(lmFitZ,2,function(tempZ){
	return(names(which(rank(-tempZ)<=nFeat)))
}))))

load("BdataBloodRefCombined.rd");#Download from https://doi.org/10.5281/zenodo.20816375
dataB<-dataB[intersect(rownames(dataB),topFeatures),]
n<-ncol(dataB)
Age<-as.numeric(phenData[,"Age"])/100

load("DNAmeEstPropV2.rd");#Output from DNAmeEstPropV2.R
nBoot<-min(nBoot,length(estPropSamps))

lmFitSamps<-lapply(1:nrow(dataB),function(tempI){
	
	tempTime<-system.time({
		
		beta<-dataB[tempI,]
				
		tempSamps<-do.call(rbind,mclapply(1:nBoot,function(iBoot){		
			tempProp<-estPropSamps[[iBoot]]
			tempWeights<-n*as.vector(rdirichlet(1,rep(1,n)))
			B<-tempProp[,"B"]
			NK<-tempProp[,"NK"]
			CD4T<-tempProp[,"CD4T"]
			CD8T<-tempProp[,"CD8T"]
			Mono<-tempProp[,"Mono"]
			Gran<-tempProp[,"Gran"]				
			tempFit<-summary(lm(as.formula("beta~0+B+NK+CD4T+CD8T+Mono+Gran+Age:B+Age:NK+Age:CD4T+Age:CD8T+Age:Mono+Age:Gran"),weights=tempWeights))[["coefficients"]][c("B:Age","NK:Age","CD4T:Age","CD8T:Age","Mono:Age","Gran:Age"),c("t value","Pr(>|t|)")]		
			tempZ<-as.vector((-qnorm(tempFit[,"Pr(>|t|)"]/2))*(2*(tempFit[,"t value"]>0)-1))
			return(tempZ)
		},mc.cores=numCores))
		colnames(tempSamps)<-c("B","NK","CD4T","CD8T","Mono","Gran")
		tempSamps[tempSamps==Inf]<-1.1*max(tempSamps[tempSamps!=Inf])
		tempSamps[tempSamps==-Inf]<-1.1*min(lmFitZ[tempSamps!=-Inf])
		
	})[[3]]	
	print(paste0(tempI,"/",nrow(dataB),"; t=",signif(tempTime,digits=2),"s"))

	return(tempSamps)

})

names(lmFitSamps)<-rownames(dataB)

save(list=c("lmFitSamps","estProp","phenData"),file="DNAmeLMageV2.rd")
