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
nFeat<-1000
libDir<-.libPaths()

CRANpackages<-c("abind")
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

load("DNAmeLMageV1.rd");#Output from DNAmeLMageV1.R
topFeatures<-apply(lmFitZ,2,function(tempZ){
	tempZ<-sort(tempZ,decreasing=T)
	return(names(head(tempZ,nFeat)))
})

load("DNAmeLMageV2.rd");#Output from DNAmeLMageV2.R
lmFitSamps<-aperm(abind(lmFitSamps,along=3),c(3,2,1))

cellTypes<-c("Gran","B","NK","CD8T","CD4T")
plotCols<-c("violetred","tan4","paleturquoise3","skyblue3","aquamarine3")
names(plotCols)<-cellTypes

plotData<-aperm(abind(sapply(cellTypes,function(tempType){
	return(sapply(1:nFeat,function(tempI){
		tempSamps<-lmFitSamps[topFeatures[tempI,tempType],tempType,]
		return(c(quantile(tempSamps,0.25,names=F),mean(tempSamps),quantile(tempSamps,0.75,names=F)))	
	}))
},simplify=F,USE.NAMES=T),along=3),c(3,1,2))
dimnames(plotData)<-list(cellTypes,c("L95","mean","U95"))


load("DNAmeLMageV3.rd");#Output from DNAmeLMageV3.R

plotData2<-aperm(abind(sapply(cellTypes,function(tempType){
	tempPlotData<-t(apply(lmFitSamps[[tempType]],2,function(tempSamps){
		return(c(quantile(tempSamps,0.25,names=F),mean(tempSamps),quantile(tempSamps,0.75,names=F)))
	}))
	colnames(tempPlotData)<-c("L95","mean","U95")
	return(tempPlotData)
},simplify=F,USE.NAMES=T),along=3),c(3,2,1))

pdf(file="DNAmeLMageCI.pdf",width=10,height=4,useDingbats=F);{
	par(mfrow=c(1,2),mar=c(3,4,1.5,1.5))
		plot(NA,xlim=c(1,1000),ylim=c(min(0,min(plotData[,"L95",])),max(plotData[,"U95",])),xlab="",ylab="",main="",xaxs="i",yaxs="i",xaxt="n",cex.axis=1.2)
	lines(100*c(1,1),c(min(0,min(plotData[,"L95",])),max(plotData[,"U95",])),lty=1,lwd=1.5)
	axis(1,at=c(1,100,1000),cex.axis=1.2)
	mtext("CpG (feature) rank",1,line=1.5,cex=1.4)
	mtext(bquote("Cell-age biomarker CpG "*italic(z)*"-score"),2,line=2.5,cex=1.4)	
	invisible(sapply(cellTypes,function(tempType){
		lines(1:nFeat,plotData[tempType,"U95",],col=plotCols[tempType],lty=2,lwd=0.7)
		lines(1:nFeat,plotData[tempType,"L95",],col=plotCols[tempType],lty=2,lwd=0.7)
		lines(1:nFeat,plotData[tempType,"mean",],col=plotCols[tempType],lty=1,lwd=1)	
		return(NULL)	
	}))

legend("topright",legend=c(names(plotCols),"Quartiles"),col=c(plotCols,"black"),lty=c(rep(1,length(plotCols)),2),lwd=2,cex=1.1)
		plot(NA,xlim=c(1,1000),ylim=c(min(0,min(plotData2[,"L95",])),max(plotData2[,"U95",])),xlab="",ylab="",main="",xaxs="i",yaxs="i",xaxt="n",cex.axis=1.2)
	lines(100*c(1,1),c(min(0,min(plotData[,"L95",])),max(plotData[,"U95",])),lty=1,lwd=1.5)
	axis(1,at=c(1,100,1000),cex.axis=1.2)
	mtext("CpG (feature) rank",1,line=1.5,cex=1.4)
	mtext(bquote("Cell-age estimator "*italic(z)*"-score"),2,line=2.5,cex=1.4)	
	invisible(sapply(cellTypes,function(tempType){
		lines(1:nFeat,plotData2[tempType,"U95",],col=plotCols[tempType],lty=2,lwd=1.5)
		lines(1:nFeat,plotData2[tempType,"L95",],col=plotCols[tempType],lty=2,lwd=1.5)
		lines(1:nFeat,plotData2[tempType,"mean",],col=plotCols[tempType],lty=1,lwd=2.5)	
		return(NULL)	
	}))

}
dev.off()
