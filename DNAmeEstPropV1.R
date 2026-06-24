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

CRANpackages<-c("MASS")
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

RLMPIfunc<-function(tempX,tempA){;#X=AW'+e
	adj<-1e-15
	commonFeat<-intersect(rownames(tempX),rownames(tempA))
	tempX<-tempX[commonFeat,]
	tempA<-tempA[commonFeat,]
	tempW<-apply(tempX,2,function(tempData){
		tempProps<-as.vector(summary(rlm(tempData~tempA,maxit=100))[["coefficients"]][-1,1])
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

estProp<-t(RLMPIfunc(dataB,centDHSbloodDMC.m));#estProp is W in X=AW'+e
estProp<-cbind(estProp,Gran=estProp[,"Neutro"]+estProp[,"Eosino"])
estProp<-estProp[,!(colnames(estProp)%in%c("Neutro","Eosino")),drop=F]

save(list=c("estProp","phenData"),file="DNAmeEstPropV1.rd")
