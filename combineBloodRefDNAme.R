rm(list=ls())


###settings and packages

dataSets<-c("GSE36054","GSE42861","GSE51032","GSE52588","GSE55763","GSE64495","GSE69270","GSE72680","GSE72775","GSE73103","GSE84727","GSE87571","GSE87648","GSE97362","GSE125367","GSE147221","EGAS00001005626")

libDir<-.libPaths()


CRANpackages<-c("MASS","rARPACK","R.utils","BiocManager")
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

BIOCpackages<-c("memoise","genefilter","sva")
for(iP in 1:length(BIOCpackages)){
	tempPackage<-BIOCpackages[iP]
	pckg0 = try(require(tempPackage,character.only=T))
	if(!pckg0){
		pckg = try(require(tempPackage,lib.loc=libDir,character.only=T))
		if(!pckg){
			print(paste0("Installing '",tempPackage,"' from bioconductor"))
			BiocManager::install(tempPackage,lib=libDir,destdir=libDir)
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


###load data

load("centDHSbloodDMC.rd");#Download from https://doi.org/10.5281/zenodo.20816375

allCpGlist<-sapply(dataSets,function(tempSet){
	print(tempSet)
	load(paste0("Bdata",tempSet,".rd"))
	if(tempSet=="GSE154566"){
		tempCpGs<-intersect(rownames(dataB450K),rownames(dataBEPIC))
	}else{
		tempCpGs<-rownames(dataB)
	}	
	return(tempCpGs)
},simplify=F,USE.NAMES=T)

allCpGs<-names(which(table(unlist(allCpGlist))==length(dataSets)))
rm(list="allCpGlist")

load("BdataGSE52588.rd")
dataB<-dataB[allCpGs,]
phenData<-phenData[colnames(dataB),]
rownames(phenData)<-phenData[,"GSM"]
colnames(dataB)<-rownames(phenData)
phenData<-cbind(Batch=rep("GSE52588",nrow(phenData)),Sample=as.character(rownames(phenData)),Age=as.character(phenData[,"Age"]))
rownames(phenData)<-phenData[,"Sample"]

dataBall<-dataB
phenDataAll<-phenData
rm(list=c("dataB","phenData"))

load("BdataEGAS00001005626.rd")
dataB<-dataB[allCpGs,]
phenData<-cbind(Batch=rep("EGAS00001005626",nrow(phenData)),Sample=as.character(rownames(phenData)),Age=as.character(phenData[,"Age_at_Consent"]))
rownames(phenData)<-phenData[,"Sample"]

dataBall<-cbind(dataBall,dataB)
phenDataAll<-rbind(phenDataAll,phenData)
rm(list=c("dataB","phenData"))

for(tempSet in setdiff(dataSets,c("GSE52588","EGAS00001005626"))){
	
	print(tempSet)

	load(paste0("Bdata",tempSet,".rd"))
	dataB<-dataB[allCpGs,]
	phenData<-cbind(Batch=rep(tempSet,nrow(phenData)),Sample=as.character(rownames(phenData)),Age=as.character(phenData[,"Age"]))
	rownames(phenData)<-phenData[,"Sample"]

	dataBall<-cbind(dataBall,dataB)
	phenDataAll<-rbind(phenDataAll,phenData)
	rm(list=c("dataB","phenData"))		

}

phenData<-phenDataAll[colnames(dataBall),]
rm(list="phenDataAll")
dataB<-dataBall
rm(list="dataBall")

dataBsvd0<-svds(dataB-rowMeans(dataB),k=100,nu=100,nv=100)
save(list=c("dataBsvd0","phenData","allCpGs"),file="BdataBloodRef0combinedSVD.rd");#Analysed by combineBloodRefDNAmePCAtest.R
rm(list="dataBsvd0")


###batch correction

dataB<-ComBat(dataB[,rownames(phenData)],batch=phenData[,"Batch"],mod=model.matrix(~Age,data=data.frame(Age=as.numeric(phenData[,"Age"]))))
dataB[dataB<0]<-0
dataB[dataB>1]<-1

estProp<-t(RLMPIfunc(dataB,centDHSbloodDMC.m))
estProp<-cbind(estProp,Gran=estProp[,"Neutro"]+estProp[,"Eosino"])
estProp<-estProp[,!(colnames(estProp)%in%c("Neutro","Eosino")),drop=F]

save(list=c("dataB","phenData","estProp"),file="BdataBloodRefCombined.rd");#Archived at https://doi.org/10.5281/zenodo.20816375

dataBsvd<-svds(dataB-rowMeans(dataB),k=100,nu=100,nv=100)
save(list=c("dataBsvd","estProp","phenData","allCpGs"),file="BdataBloodRefCombinedSVD.rd");#Analysed by combineBloodRefDNAmePCAtest.R
