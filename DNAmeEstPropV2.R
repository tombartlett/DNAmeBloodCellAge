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
libDir<-.libPaths()

CRANpackages<-c("MASS","parallel","MCMCpack")
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

RLMPIfunc<-function(tempX,tempA,tempW){;#X=AW'+e
	adj<-1e-15
	commonFeat<-intersect(rownames(tempX),rownames(tempA))
	tempX<-tempX[commonFeat,]
	tempA<-tempA[commonFeat,]
	tempW<-apply(tempX,2,function(tempData){
		tempProps<-as.vector(summary(rlm(tempData~tempA,weights=tempW,wt.method="case",maxit=100))[["coefficients"]][-1,1])
		tempProps[tempProps<0]<-0;
		tempProps<-tempProps/sum(tempProps)
		return(tempProps)
	})
	rownames(tempW)<-colnames(tempA)
	colnames(tempW)<-colnames(tempX)
	return(tempW)
}

load("centDHSbloodDMC.rd");#Download from https://doi.org/10.5281/zenodo.20816375

load("BdataEstPropBloodRefCombined.rd");#Download from https://doi.org/10.5281/zenodo.20816375
dataB<-dataBestProp
rm(list="dataBestProp")
n<-nrow(dataB)

estPropSamps<-mclapply(1:nBoot,function(iBoot){
	tempTime<-system.time({
		tempWeights<-n*as.vector(rdirichlet(1,rep(1,n)))	
		estPropTemp<-t(RLMPIfunc(dataB,centDHSbloodDMC.m,tempWeights));#estProp is W in X=AW'+e
		estPropTemp<-cbind(estPropTemp,Gran=estPropTemp[,"Neutro"]+estPropTemp[,"Eosino"])
		estPropTemp<-estPropTemp[,!(colnames(estPropTemp)%in%c("Neutro","Eosino")),drop=F]
	})[[3]]
	print(paste0(iBoot,"/",nBoot,"; t=",signif(tempTime,digits=2),"s"))
	return(estPropTemp)
},mc.cores=numCores)

save(list=c("estPropSamps","phenData"),file="DNAmeEstPropV2.rd")
