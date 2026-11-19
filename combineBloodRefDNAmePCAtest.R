rm(list=ls())

nPC<-20
libDir<-.libPaths()

CRANpackages<-c("Hmisc")
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

load("BdataBloodRef0combinedSVD.rd")
batchCor0<-abs(sapply(unique(phenData[,"Batch"]),function(tempBatch){
	inTemp<-as.vector(phenData[,"Batch"])%in%tempBatch
	return(apply(dataBsvd0[["v"]][,1:nPC,drop=F],2,function(tempVals){
		return(somers2(tempVals,inTemp)[["Dxy"]])
	}))	
}))
ageCor0<-abs(as.vector(cor(as.numeric(phenData[,"Age"]),dataBsvd0[["v"]][,1:nPC,drop=F])))

load("BdataBloodRefCombinedSVD.rd")
batchCor<-abs(sapply(unique(phenData[,"Batch"]),function(tempBatch){
	inTemp<-as.vector(phenData[,"Batch"])%in%tempBatch
	return(apply(dataBsvd[["v"]][,1:nPC,drop=F],2,function(tempVals){
		return(somers2(tempVals,inTemp)[["Dxy"]])
	}))	
}))
ageCor<-abs(as.vector(cor(as.numeric(phenData[,"Age"]),dataBsvd[["v"]][,1:nPC,drop=F])))

plotData<-rbind(Age=ageCor,t(batchCor))
colnames(plotData)<-paste(1:ncol(plotData))
plotData<-t(plotData[rev(rownames(plotData)),])
xLoc<-1:nrow(plotData)
yLoc<-1:ncol(plotData)
heatCols<-c(colorRampPalette(c("white","goldenrod3"))(6)[-6],"goldenrod3",colorRampPalette(c("goldenrod3","tomato3"))(6)[-1])
heatBreaks<-seq(0,1.1,by=0.1)-0.05
heatVals<-seq(0,1,by=0.1)
pdf(file="combineBloodRefDNAmePCAtest.pdf",width=9.5,height=7.5);
	par(mar=c(5,12.5,4,1))
	image(x=xLoc,y=yLoc,z=plotData,col=heatCols,breaks=heatBreaks,xaxt="n",yaxt="n",xlab="",ylab="",main="")
	axis(1,at=xLoc,labels=rownames(plotData),lwd=1.5,tck=-0.015,cex.axis=1.5,las=2)
	axis(2,at=yLoc,labels=colnames(plotData),lwd=1.5,tck=-0.015,cex.axis=1.5,las=1)
	mtext("PC",side=1,line=3.5,cex=1.8)
	par(xpd=NA)
	legend(grconvertX(0.5,from="npc", to="user"),grconvertY(1,from="npc",to="user"),legend=paste0(heatVals),fill=heatCols,bty="n",xjust=0.5,yjust=0,horiz=T,x.intersp=0.4,cex=1.4,title="Correlation",text.width=0.9*strwidth(paste0(heatVals)),title.cex=1.5)
	rm(list="plotData")
dev.off()

plotData<-rbind(Age=ageCor0,t(batchCor0))
colnames(plotData)<-paste(1:ncol(plotData))
plotData<-t(plotData[rev(rownames(plotData)),])
xLoc<-1:nrow(plotData)
yLoc<-1:ncol(plotData)
heatCols<-c(colorRampPalette(c("white","goldenrod3"))(6)[-6],"goldenrod3",colorRampPalette(c("goldenrod3","tomato3"))(6)[-1])
heatBreaks<-seq(0,1.1,by=0.1)-0.05
heatVals<-seq(0,1,by=0.1)
pdf(file="combineBloodRefDNAmePCAtest0.pdf",width=9.5,height=7.5);
	par(mar=c(5,12.5,4,1))
	image(x=xLoc,y=yLoc,z=plotData,col=heatCols,breaks=heatBreaks,xaxt="n",yaxt="n",xlab="",ylab="",main="")
	axis(1,at=xLoc,labels=rownames(plotData),lwd=1.5,tck=-0.015,cex.axis=1.5,las=2)
	axis(2,at=yLoc,labels=colnames(plotData),lwd=1.5,tck=-0.015,cex.axis=1.5,las=1)
	mtext("PC",side=1,line=3.5,cex=1.8)
	par(xpd=NA)
	legend(grconvertX(0.5,from="npc", to="user"),grconvertY(1,from="npc",to="user"),legend=paste0(heatVals),fill=heatCols,bty="n",xjust=0.5,yjust=0,horiz=T,x.intersp=0.4,cex=1.4,title="Correlation",text.width=0.9*strwidth(paste0(heatVals)),title.cex=1.5)
	rm(list="plotData")
dev.off()

pdf(file="combineBloodRefDNAmePCPVE.pdf",height=8,width=8);{
	par(mar=c(5,5,1,1))
	tempD2<-dataBsvd[["d"]]^2
	Lvars<-tempD2/sum(tempD2)
	tempLoc<-c(1,seq(10,length(Lvars),by=10))
	plot(Lvars,pch=20,col="black",ylim=c(0,0.5),xlab="",ylab="",main="",cex.axis=1.5,xaxt="n")
	axis(1,at=tempLoc,labels=paste(tempLoc),lwd=1.5,tck=-0.015,cex.axis=1.5,las=1)
	mtext("PC",side=1,line=3.2,cex=1.8)
	mtext("PVE",side=2,line=3.5,cex=1.8)
	rm(list=c("Lvars","tempD2","tempLoc"))
}
dev.off() 

pdf(file="combineBloodRefDNAmePCPVE0.pdf",height=8,width=8);{
	par(mar=c(5,5,1,1))
	tempD2<-dataBsvd0[["d"]]^2
	Lvars<-tempD2/sum(tempD2)
	tempLoc<-c(1,seq(10,length(Lvars),by=10))
	plot(Lvars,pch=20,col="black",ylim=c(0,0.5),xlab="",ylab="",main="",cex.axis=1.5,xaxt="n")
	axis(1,at=tempLoc,labels=paste(tempLoc),lwd=1.5,tck=-0.015,cex.axis=1.5,las=1)
	mtext("PC",side=1,line=3.2,cex=1.8)
	mtext("PVE",side=2,line=3.2,cex=1.8)
	rm(list=c("Lvars","tempD2","tempLoc"))
}
dev.off() 
