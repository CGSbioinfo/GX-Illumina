library(lumi)
library(limma)
#library(hwriter) not needed v1.01
library(ReportingTools)
library(ggplot2)
library(gplots)
library(reshape)
library(annotate)
library(illuminaMousev2.db)
setClass("IlluminaData",
         contains="LumiBatch",
         representation=representation(
           designMatrix="matrix",
           groups="data.frame",
           raw_data_file="character",
           norm_data="matrix",
           results="data.frame",
           controlData="data.frame"
           )
         )

setMethod("initialize","IlluminaData",function(.Object,...,rawDataFile="", groupFile=NA, results=NA, .volcano=NA){
  if(is.na(groupFile)){
    .design=NA
    .Object@groups<-data.frame()
  }else{
    .groups=read.table(groupFile, sep="\t", header=T)
    .Object@groups<-.groups
  }
        
  if(is.na(results)){
    .Object@results<-data.frame()
  }
  
            
  .Object@raw_data_file<-rawDataFile
  if(rawDataFile != ""){
    .lumiObject<-lumiR(rawDataFile)
    if(is.na(groupFile)){
      .lumiObject.set<-.lumiObject
    }else{
      
      .lumiObject.set<-.lumiObject[,as.character(.groups$SampleID)]
      #print(sampleNames(.lumiObject))
      .groups<-.groups[match(sampleNames(.lumiObject.set),.groups$SampleID),]
      #print(.groups)
      
    }
    .Object@assayData<-assayData(.lumiObject.set)
    .Object@phenoData<-phenoData(.lumiObject.set)
    .Object@experimentData<-experimentData(.lumiObject.set)
    .Object@annotation<-annotation(.lumiObject.set)
    .Object@controlData<-controlData(.lumiObject.set)
    .Object@protocolData<-.Object@phenoData
    pData(featureData(.Object))<-pData(featureData(.lumiObject.set))
    .Object@groups<-.groups
              
  }
  
  .Object@designMatrix<-matrix()
  return(.Object)
})

loadIllumina<-function(rawDataFile,groupFile=NA){
  object<-new("IlluminaData",rawDataFile=rawDataFile,groupFile=groupFile)
  return(object)
}

setMethod("show", "IlluminaData",
          function(object){
            print(paste("Number of features: ",dim(exprs(object))[1],sep=""))
            print(paste("Number of samples: ",dim(exprs(object))[2],sep=""))
            print(paste("Number of Groups: ",as.character(length(levels(as.factor(groups(object)[,2])))),sep=""))
            print("Group Names: ")
            print(as.character(levels(as.factor(groups(object)[,2]))))
          })

setGeneric("technicalRplAve",function(object,group=NA) standardGeneric("technicalRplAve"))
setMethod("technicalRplAve", "IlluminaData",
          function(object,group){
            group.info=cbind(group,as.factor(paste(group$Group,group$Sibship,sep="_")))
            colnames(group.info)[ncol(group.info)]="block"
            block = as.character(unique(group.info$block))
            m=as.data.frame(matrix(nrow= nrow(exprs(object)),ncol=length(block)))
            rownames(m)=rownames(exprs(object))
            m.se=as.data.frame(matrix(nrow= nrow(exprs(object)),ncol=length(block)))
            rownames(m.se)=rownames(exprs(object))
            m.detect=as.data.frame(matrix(nrow= nrow(exprs(object)),ncol=length(block)))
            rownames(m.detect)=rownames(exprs(object))
            for (i in 1:length(block)){
              samples=as.character(group.info$SampleID[group.info$block==block[i]])
              # Mean intensity values
              data=exprs(object[,samples])
              norm.data=normalize.quantiles(data)
              rownames(norm.data)=rownames(data)
              data=rowMeans(norm.data)
              m[,i]=data[rownames(m)]
              colnames(m)[i]=unique(gsub("_[0-9]","",samples)) ## sample names of replicates should end as _1, _2, _3, ...
              
              # Mean se values
              data=se.exprs(object[,samples])
              data=rowMeans(data)
              m.se[,i]=data[rownames(m.se)]
              colnames(m.se)[i]=unique(gsub("_[0-9]","",samples)) ## sample names of replicates should end as _1, _2, _3, ...
              
              #Mean detection pval values
              data=detection(object[,samples])
              data=rowMeans(data)
              m.detect[,i]=data[rownames(m.detect)]
              colnames(m.detect)[i]=unique(gsub("_[0-9]","",samples)) ## sample names of replicates should end as _1, _2, _3, ...
            }
            exprs(object)=as.matrix(m)
            se.exprs(object)=as.matrix(m.se)
            detection(object)=as.matrix(m.detect)
            new.groups=as.data.frame(cbind(SampleID=colnames(m), Group=c(unique(as.character(group.info$Group))), 
                                           Sibship=rep(c(unique(group.info$Sibship)),each = 2)))
            object@groups=new.groups
            object@groups<-droplevels(object@groups)
            phenoData(object)@data=as.data.frame(matrix(as.character(colnames(m))))
            rownames(phenoData(object)@data)=colnames(m)
            colnames(phenoData(object)@data)="sampleID"
            return(object)
          })

setGeneric("processData",function(object,group=NA, rpl.mean=FALSE) standardGeneric("processData"))
setMethod("processData", "IlluminaData",
          function(object,group,rpl.mean){
            if(all(is.na(group))){
              object<-object[detectionCall(object)>0,]
              object<-lumiT(object)
              object<-lumiN(object)
            }else{
              object<-object[,as.character(group[,1])]
              object<-object[detectionCall(object)>0,]
              if (rpl.mean==TRUE){
                object=technicalRplAve(object,group)
                print("rpl.mean")
              }
              object<-lumiT(object)
              object<-lumiN(object)
              if (!rpl.mean==TRUE){
                object@groups<-group
                #for subsetting it is needed to drop levels that are not used, otherwise they will interfere with 
                #desing matrix creation
                object@groups<-droplevels(object@groups)
              }  
            }
            return(object)
          })

setGeneric("addDetectToNormData",function(object) standardGeneric("addDetectToNormData"))
setMethod("addDetectToNormData", "IlluminaData",
          function(object){
            normData.exprs=as.data.frame(exprs(object))
            normData.exprs=cbind(normData.exprs,detectionCall(object, type="probe"))
            row.order=match(row.names(normData.exprs), row.names(detection(object)))
            detect_table=as.data.frame(detection(object)[row.order,])
            names=colnames(detect_table)
            for (k in 1:length(names)){
              index<-which(colnames(normData.exprs)==names[k])
              name=names[k]
              normData.exprs=cbind(normData.exprs[,1:index], detect_table[,names[k]]<=0.01,
                                  detect_table[,names[k]], normData.exprs[,(index+1):dim(normData.exprs)[2]])
              colnames(normData.exprs)[((index):(index+2))]=c(names[k],paste(names[k],".Detected",sep=""),paste(names[k],".Detect.pval",sep=""))
            }
            colnames(normData.exprs)[(dim(normData.exprs)[2])]="Overall.detection"
            return(normData.exprs)
          })


setGeneric("setAnnotation", function(object,arrayType) standardGeneric("setAnnotation"))
setMethod("setAnnotation", "data.frame",
          function(object,arrayType=""){
            if(arrayType=="Mouse"){
              if(require("illuminaMousev2.db")){
                idList=as.character(rownames(object))
                mappedID<-as.character(lookUp(idList, "illuminaMousev2.db", "SYMBOL"))
                temp<-data.frame(ID=idList,SYMBOL=mappedID)
                temp[temp=="NA"] <- NA
                temp$ID<-as.character(temp$ID)
                order.ID=match(rownames(object),temp[,1])
                temp<-as.data.frame(temp[order.ID,])
                temp<-cbind(object,temp$SYMBOL)
                colnames(temp)[dim(temp)[2]]<-"GENE_SYMBOL"
              }
            } else if(arrayType=="Human"){
              if(require("illuminaHumanv4.db")){
                idList=as.character(rownames(object))
                x <- illuminaHumanv4SYMBOLREANNOTATED
                xx <- unlist(as.list(x[idList]))
                index=match(idList,names(xx))
                temp=cbind(object,xx[index])
                colnames(temp)[ncol(temp)]="Gene_Symbol"
                #mappedID<-as.character(lookUp(idList, "illuminaHumanv4.db", "SYMBOL"))
                #temp<-data.frame(ID=idList,SYMBOL=mappedID)
                #temp[temp=="NA"] <- NA
                #temp$ID<-as.character(temp$ID)
                #order.ID=match(rownames(object),temp[,1])
                #temp<-as.data.frame(temp[order.ID,])
                #temp<-cbind(object,temp$SYMBOL)
                #colnames(temp)[dim(temp)[2]]<-"GENE_SYMBOL"
                
                #geneID
                x <- illuminaHumanv4ENTREZREANNOTATED
                xx <- unlist(as.list(x[idList]))
                index=match(rownames(temp),names(xx))
                temp=cbind(temp,xx[index])
                colnames(temp)[ncol(temp)]="Entrez_ID"
                
                # Chr
                x <- illuminaHumanv4CHR
                xx <- unlist(as.list(x[idList]))
                index=match(rownames(temp),names(xx))
                temp=cbind(temp,xx[index])
                colnames(temp)[ncol(temp)]="chr"
                
                # GENOMIC LOCATION
                x <- illuminaHumanv4GENOMICLOCATION
                # Convert to a list
                xx <- unlist(as.list(x[idList]))
                index=match(rownames(temp),names(xx))
                temp=cbind(temp,xx[index])
                colnames(temp)[ncol(temp)]="Probe_Genomic_Coordinates_hg19"
                
                # PROBE SEQUENCE
                x <- illuminaHumanv4PROBESEQUENCE
                # Convert to a list
                xx <- unlist(as.list(x[idList]))
                index=match(rownames(temp),names(xx))
                temp=cbind(temp,xx[index])
                colnames(temp)[ncol(temp)]="Probe_sequence"
                
                # PROBE QUALITY
                x <- illuminaHumanv4PROBEQUALITY
                # Convert to a list
                xx <- unlist(as.list(x[idList]))
                index=match(rownames(temp),names(xx))
                temp=cbind(temp,xx[index])
                colnames(temp)[ncol(temp)]="Probe quality"
                
                # PROBE Similarity to genomic match
                #x <- illuminaHumanv4GENOMICMATCHSIMILARITY
                ## Convert to a list
                #xx <- unlist(as.list(x[idList]))
                #index=match(rownames(temp),names(xx))
                #temp=cbind(temp,xx[index])
                #colnames(temp)[ncol(temp)]="Probe_match_similarity"
                
              }} else{
                print("Error, not arrayType specified, requiere Mouse or Human") 
              }
              return(temp)
            })


### generatedHybControlGraph added v1.01
setGeneric("generateHybControlGraph", function(object, path,arrayType) standarGeneric("generateHybControlGraph"))
setMethod("generateHybControlGraph","IlluminaData",function(object,path="character",arrayType="character"){
  l=0
  h=0
  k=0
  controlTable<-getControlData(object)
  cy3_control<-controlTable[controlTable[,1]=="CY3_HYB",]
  if(arrayType=="human"){
    for(i in 1:6){
      if(cy3_control[i,]$ProbeID=="4010327" || cy3_control[i,]$ProbeID=="1110170"){
        l=l+1
        cy3_control[i,]$controlType<-paste("High_",l,sep="")
      }
      else if(cy3_control[i,]$ProbeID=="2510500" || cy3_control[i,]$ProbeID=="4610291"){
        h=h+1
        cy3_control[i,]$controlType<-paste("Medium_",h,sep="")
      }
      else if(cy3_control[i,]$ProbeID=="7560739" || cy3_control[i,]$ProbeID=="1450438"){
        k=k+1
        cy3_control[i,]$controlType<-paste("Low_",k,sep="")
      }
    }
  }else if(arrayType=="mouse"){
    for(i in 1:6){
      if(cy3_control[i,]$ProbeID=="7200743" || cy3_control[i,]$ProbeID=="6450180" || cy3_control[i,]$ProbeID=="3800196"){
        l=l+1
        cy3_control[i,]$controlType<-paste("High_",l,sep="")
      }
      else if(cy3_control[i,]$ProbeID=="1400044" || cy3_control[i,]$ProbeID=="2600040"){
        h=h+1
        cy3_control[i,]$controlType<-paste("Medium_",h,sep="")
      }
      else if(cy3_control[i,]$ProbeID=="730475" || cy3_control[i,]$ProbeID=="5820544"){
        k=k+1
        cy3_control[i,]$controlType<-paste("Low_",k,sep="")
      }
    }
  }
  melt<-melt(data=cy3_control[,c(1,3:length(controlTable))], id.vars="controlType", variable_name="Samples")
  
  #plot using ggplot
  hybPlot<-ggplot(melt, aes(x=Samples ,y=value,colour=controlType),las=2)+geom_point(aes(group=Samples))+theme(axis.text.x=element_text(angle=90))+ylab("Intensity values")  
  png(filename=paste(path,"/Hyb_Controls.png", sep=""), width = 700, height = 700)
  print(hybPlot)
  graphics.off()
})

#generateStringencyControlGraph added v1.01
setGeneric("generateStringencyControlGraph", function(object, path, arrayType) standarGeneric("generateStringencyControlGraph"))
setMethod("generateStringencyControlGraph","IlluminaData",function(object,path="character",arrayType="character"){

  h=0
  l=0
  sum_pm=0
  sum_mm=0
  controlTable<-getControlData(object)
  low_string_control<-controlTable[controlTable[,1]=="LOW_STRINGENCY_HYB",]
  #renaming the data according to probe id
  if (arrayType=="human"){
    for(i in 1:8){
      if(low_string_control[i,]$ProbeID=="1110170" || low_string_control[i,]$ProbeID=="2510500" || low_string_control[i,]$ProbeID=="4010327" || low_string_control[i,]$ProbeID=="4610291"){
        l=l+1
        low_string_control[i,]$controlType<-paste("PM_",l,sep="")
        sum_pm=sum_pm+low_string_control[i,3:length(low_string_control[i,])]
      }
      else{
        h=h+1
        low_string_control[i,]$controlType<-paste("MM_",h,sep="")
        sum_mm=sum_mm+low_string_control[i,3:length(low_string_control[i,])]
      }
    }
  }else if(arrayType=="mouse"){
    for(i in 1:8){
        if(low_string_control[i,]$ProbeID=="2600040" || low_string_control[i,]$ProbeID=="3800196" || low_string_control[i,]$ProbeID=="1400044" || low_string_control[i,]$ProbeID=="6450180"){
            l=l+1
            low_string_control[i,]$controlType<-paste("PM_",l,sep="")
            sum_pm=sum_pm+low_string_control[i,3:length(low_string_control[i,])]
        }else{
            h=h+1
            low_string_control[i,]$controlType<-paste("MM_",h,sep="")
            sum_mm=sum_mm+low_string_control[i,3:length(low_string_control[i,])]
        }
      }
  }
  mean_pm=cbind("Perfect match","mean",sum_pm/4)
  colnames(mean_pm)[1:2]=c("controlType", "ProbeID")
  mean_mm=cbind("Mismatch","mean",sum_mm/4)
  colnames(mean_mm)[1:2]=c("controlType", "ProbeID")
  low_string_control=rbind(low_string_control,mean_pm,mean_mm)
  melt<-melt(data=low_string_control[c(9:10),c(1,3:length(controlTable))], id.vars="controlType", variable_name="Samples") # Plot mean value  

  #plot using ggplot
  stingPlot<-ggplot(melt, aes(x=Samples ,y=value,colour=controlType),las=2)+geom_point(aes(group=Samples))+theme(axis.text.x=element_text(angle=90))+ylab("Intensity values")
  png(filename=paste(path,"/String_Controls.png",sep=""), width = 700, height = 700)
  print(stingPlot)
  graphics.off()
})


#generateLabelControlGraph added v1.01
setGeneric("generateLabelControlGraph", function(object, path) standarGeneric("generateLabelControlGraph"))
setMethod("generateLabelControlGraph","IlluminaData",function(object,path="character"){
  controlTable<-getControlData(object)
  label_control<-controlTable[controlTable[,1]=="BIOTIN",]
  negative<-controlTable[controlTable[,1]=="NEGATIVE",]
  merged<-rbind(label_control, negative[1:10,])
  
  #melting the data
  for(i in 1:12){
    newName<-paste(merged[i,1],"_",i, sep="")
    merged[i,1]<-newName
  }
  melt<-melt(data=merged[,c(1,3:length(controlTable))], id.vars="controlType", variable_name="Samples")
  #plot using ggplot
  labelPlot<-ggplot(melt, aes(x=Samples ,y=value,colour=controlType),las=2)+geom_point(aes(group=Samples))+theme(axis.text.x=element_text(angle=90))+ylab("Intensity values")
  png(filename=paste(path,"/Label_Controls.png",sep=""), width = 700, height = 700)
  print(labelPlot)
  graphics.off()
})

#generateNegativeControlGraph added v1.01
setGeneric("generateNegativeControlGraph", function(object, path) standarGeneric("generateNegativeControlGraph"))
setMethod("generateNegativeControlGraph","IlluminaData",function(object,path="character"){
  controlTable<-getControlData(object) 
  merged<-controlTable[controlTable[,1]=="NEGATIVE",]
   
  #melting the data
  for(i in 1:12){
    newName<-paste(merged[i,1],"_",i, sep="")
    merged[i,1]<-newName
  }
  melt<-melt(data=merged[,c(1,3:length(controlTable))], id.vars="controlType", variable_name="Samples")
  #plot using ggplot
  negPlot<-ggplot(melt, aes(x=Samples ,y=value,colour=controlType),las=2)+geom_point(aes(group=Samples))+theme(axis.text.x=element_text(angle=90))+ylab("Intensity values")
  png(filename=paste(path,"/Neg_Controls.png",sep=""), width = 700, height = 700)
  print(negPlot)
  graphics.off()
  
})



## Generate Technical QC graphs added v1.01
setGeneric("generateTechQCgraph", function(object, arrayType) standardGeneric("generateTechQCgraph"))
setMethod("generateTechQCgraph","IlluminaData",function(object, arrayType="character"){
  baseDir=getwd()
  folder="Reports/QCReports"
  folderGraph2="QCGraphs"
  folderGraph=paste(getwd(),"/Reports/QCReports/QCGraphs",sep="")
  newDir=file.path(getwd(),folder)
  if(!(folder %in% dir())){
    dir.create(newDir,recursive = T)
  }
  if(!(folderGraph %in% dir(newDir))){
    dir.create(folderGraph,recursive = T)
  }
  generateHybControlGraph(object, path=folderGraph, arrayType=arrayType)
  generateStringencyControlGraph(object,path=folderGraph,arrayType=arrayType)
  generateNegativeControlGraph(object, path=folderGraph)
  generateLabelControlGraph(object,path=folderGraph)
  png(filename=paste(folderGraph,"/HC_all_samples.png",sep=""), width = 700, height = 700)
  plotSampleRelation(object)
  dev.off()
  png(filename=paste(folderGraph,"/MDS_samples.png",sep=""), width = 700, height = 700)
  plotSampleRelation(object, method="mds",color=groups(object)$Group)
  dev.off()
  png(filename=paste(folderGraph,"/Outlier_detection.png",sep=""), width = 700, height = 700)
  detectOutlier(object, ifPlot=T)
  dev.off()
  png(filename=paste(folderGraph,"/Boxplot.png",sep=""), width = 700, height = 700)
  plot(object, what="b")
  dev.off()
  png(filename=paste(folderGraph,"/Density.png",sep=""), width = 700, height = 700)
  plot(object, what="d")
  dev.off()
  
 
  
  
  ############################# Mouse array not available anymore v1.01
#   else if(arrayType=="mouse"){
#     for(i in 1:6){
#       if(cy3_control[i,]$ProbeID=="7200743" || cy3_control[i,]$ProbeID=="6450180" || cy3_control[i,]$ProbeID=="3800196"){
#         l=l+1
#         cy3_control[i,]$controlType<-paste("High_",l,sep="")
#       }
#       else if(cy3_control[i,]$ProbeID=="1400044" || cy3_control[i,]$ProbeID=="2600040"){
#         h=h+1
#         cy3_control[i,]$controlType<-paste("Medium_",h,sep="")
#       }
#       else if(cy3_control[i,]$ProbeID=="730475" || cy3_control[i,]$ProbeID=="5820544"){
#         k=k+1
#         cy3_control[i,]$controlType<-paste("Low_",k,sep="")
#       }
#     }
#     melt<-melt(data=cy3_control[,c(1,3:length(controlTable))], id.vars="controlType", variable_name="Samples")
#     #plot using ggplot
#     hybPlot<-ggplot(melt, aes(x=Samples ,y=value,colour=controlType),las=2)+geom_point(aes(group=Samples))+theme(axis.text.x=element_text(angle=90))+ylab("Intensity values")  
#     png(filename=paste(folderGraph,"/Hyb_Controls.png", sep=""), width = 700, height = 700)
#     print(hybPlot)
#     dev.off()
#     
#     
#     ##################
#     #Strigency control
#     low_string_control<-controlTable[controlTable[,1]=="LOW_STRINGENCY_HYB",]
#     
#     #melting the data
#     h=0
#     l=0
#     sum_pm=0
#     sum_mm=0
#     #renaming the data according to probe id
#     for(i in 1:8){
#       if(low_string_control[i,]$ProbeID=="2600040" || low_string_control[i,]$ProbeID=="3800196" || low_string_control[i,]$ProbeID=="1400044" || low_string_control[i,]$ProbeID=="6450180"){
#         l=l+1
#         low_string_control[i,]$controlType<-paste("PM_",l,sep="")
#         sum_pm=sum_pm+low_string_control[i,3:length(low_string_control[i,])]
#       }
#       else{
#         h=h+1
#         low_string_control[i,]$controlType<-paste("MM_",h,sep="")
#         sum_mm=sum_mm+low_string_control[i,3:length(low_string_control[i,])]
#       }
#     }
#     mean_pm=cbind("Perfect match","mean",sum_pm/4)
#     colnames(mean_pm)[1:2]=c("controlType", "ProbeID")
#     mean_mm=cbind("Mismatch","mean",sum_mm/4)
#     colnames(mean_mm)[1:2]=c("controlType", "ProbeID")
#     low_string_control=rbind(low_string_control,mean_pm,mean_mm)
#     melt<-melt(data=low_string_control[c(9:10),c(1,3:length(controlTable))], id.vars="controlType", variable_name="Samples") # Plot mean value  
#     # melt<-melt(data=low_string_control[,c(1,3:length(controlTable))], id.vars="controlType", variable_name="Samples")  Plot values separately
#     # plot using ggplot
#     stingPlot<-ggplot(melt, aes(x=Samples ,y=value,colour=controlType),las=2)+geom_point(aes(group=Samples))+theme(axis.text.x=element_text(angle=90))+ylab("Intensity values")
#     png(filename=paste(folderGraph,"/String_Controls.png",sep=""), width = 700, height = 700)
#     print(stingPlot)
#     dev.off()
#     
#     ##############
#     #labelling/neg
#     label_control<-controlTable[controlTable[,1]=="BIOTIN",]
#     negative<-controlTable[controlTable[,1]=="NEGATIVE",]
#     merged<-rbind(label_control, negative[1:10,])
#     
#     #melting the data
#     for(i in 1:12){
#       newName<-paste(merged[i,1],"_",i, sep="")
#       merged[i,1]<-newName
#     }
#     melt<-melt(data=merged[,c(1,3:length(controlTable))], id.vars="controlType", variable_name="Samples")
#     #plot using ggplot
#     labelPlot<-ggplot(melt, aes(x=Samples ,y=value,colour=controlType),las=2)+geom_point(aes(group=Samples))+theme(axis.text.x=element_text(angle=90))+ylab("Intensity values")
#     png(filename=paste(folderGraph,"/Label_Controls.png",sep=""), width = 700, height = 700)
#     print(labelPlot)
#     dev.off()
#     
#   }

###commented out temp 1.01
#   png(filename=paste(folderGraph,"/HC_all_samples.png",sep=""), width = 700, height = 700)
#   plotSampleRelation(object)
#   dev.off()
#   png(filename=paste(folderGraph,"/MDS_samples.png",sep=""), width = 700, height = 700)
#   plotSampleRelation(object, method="mds",color=groups(object)$Group)
#   dev.off()
#   png(filename=paste(folderGraph,"/Outlier_detection.png",sep=""), width = 700, height = 700)
#   detectOutlier(object, ifPlot=T)
#   dev.off()
#   png(filename=paste(folderGraph,"/Boxplot.png",sep=""), width = 700, height = 700)
#   plot(object, what="b")
#   dev.off()
#   png(filename=paste(folderGraph,"/Density.png",sep=""), width = 700, height = 700)
#   plot(object, what="d")
#   dev.off()
})

### Generate group graphs for QC added v1.01
setGeneric("generateGroupGraphs", function(object,reportName) standardGeneric("generateGroupGraphs"))
setMethod("generateGroupGraphs", "IlluminaData", function(object,reportName="character"){
  
  baseDir=getwd()
  folder=paste("Reports/",reportName,sep="")
  folderGraph2="GroupGraphs"
  folderGraph=paste(getwd(),"/Reports/",reportName,"/GroupGraphs",sep="")
  newDir=file.path(getwd(),folder)
  if(!(folder %in% dir())){
    dir.create(newDir)
  }
  if(!(folderGraph %in% dir(newDir))){
    dir.create(folderGraph)
  }
  
  group_cat<-levels(factor(groups(object)$Group))
  groups_all<-groups(object)
  if( length(group_cat)==1){
    n=1
  } else {
    n=(length(group_cat)-1)
  }
  for(i in 1:n){
    for(l in (i+1):length(group_cat)){
      if (all(is.na(groups_all[group_cat[l]==groups_all$Group,]))){
        group_selected<-groups_all
      } else {
        group_selected<-rbind(groups_all[group_cat[i]==groups_all$Group,], groups_all[group_cat[l]==groups_all$Group,]) 
      }
      
      #scatter plots
      nameFile=paste(folderGraph,"/scatter_plot_groups_",group_cat[i],"_",group_cat[l],".png",sep="")
      nameOnly=paste("scatter_plot_groups_",group_cat[i],"_",group_cat[l],".png",sep="")
      png(filename = nameFile,width=1024,height=748)
      plot(object[,as.character(group_selected$SampleID)], what="p", subset=NULL)
      dev.off()
      
      #boxplot
      nameFile=paste(folderGraph,"/boxplot_groups_",group_cat[i],"_",group_cat[l],".png",sep="")
      nameOnly=paste("boxplot_groups_",group_cat[i],"_",group_cat[l],".png",sep="")
      png(filename = nameFile)
      plot(object[,as.character(group_selected$SampleID)], what="b", subset=NULL)
      dev.off()
      
      #hierarchical clustering
      nameFile=paste(folderGraph,"/HC_groups_",group_cat[i],"_",group_cat[l],".png",sep="")
      nameOnly=paste("HC_groups_",group_cat[i],"_",group_cat[l],".png",sep="")
      png(filename = nameFile)
      plotSampleRelation(object[,as.character(group_selected$SampleID)])
      dev.off()
      if (n==1){
        break
      }
    } 
  }
})

### Generate Comparison and comparisons graphs added v1.01
setGeneric("multipleComparisons", function(object, arrayType, rpl.mean=FALSE, pairedTest=FALSE, heatmap.geneList=NA) standardGeneric("multipleComparisons"))
setMethod("multipleComparisons", "IlluminaData",
          function(object, arrayType="", rpl.mean, pairedTest, heatmap.geneList=NA){
            if(is.na(groups(object))){
              print("Groups is empty, needed to do any comparisons")
            }else{
              #png()
              results_sum<-data.frame(matrix(ncol=4))
              colnames(results_sum)<-c("Comparisons","Number of selected probe","Significant results at 0.01", "Significant results at 0.05")
              listVolca<-vector(mode="list")
              listHeatMap<-vector(mode="list")
              compName<-vector(mode="list")
              fileNameList<-vector(mode="list")
              fileList<-vector(mode="list")
              fileNameList_hm<-vector(mode="list")
              fileList_hm<-vector(mode="list")
              levels_groups=levels(as.factor(groups(object)[,2]))
              group_all=groups(object)
              for(i in 1:(length(levels_groups)-1)){
                for(l in (i+1):length(levels_groups)){
                  #create folder
                  name1=levels_groups[i]
                  name2=levels_groups[l]
                  newDir=file.path(getwd(),paste("Comparisons/",name1,"_vs_",name2,sep=""))
                  #print(newDir)
                  dir.create(newDir,recursive=T)
                  dir.create(paste(newDir,"/Graphs", sep=""))
                  
                  subgroup<-group_all[group_all[,2]==levels_groups[i] | group_all[,2]==levels_groups[l],]
                  if (rpl.mean==TRUE){
                    subdata<-processData(object,subgroup, rpl.mean=TRUE)
                    print("rpl")
                  } else {
                    subdata<-processData(object,subgroup)
                  }
                  if (pairedTest==TRUE){
                    subdata.m<-createDesign(subdata,pairedTest = TRUE)
                    subdata.comp<-compareGroup(subdata.m, pairedTest=TRUE)
                    print("paired")
                  } else {
                    subdata.m<-createDesign(subdata) 
                    subdata.comp<-compareGroup(subdata.m)
                  }
                  results_table<-as.data.frame(results(subdata.comp))
                  results_table<-results_table[order(results_table$adj.P.Val),]
                  colnames(results_table)[10]="logFC"
                  colnames(results_table)[14]="adj.P.Val"
                  
                  
                  temp=match(row.names(results_table),row.names(detection(subdata)))
                  detect_table=as.data.frame(detection(subdata)[temp,])
                  names=colnames(detect_table)
                  #Overall.detection=matrix(0, nrow = nrow(detect_table))
                  for (k in 1:length(names)){
                    index<-which(colnames(results_table)==names[k])
                    results_table=cbind(results_table[,1:index], detect_table[,names[k]]<=0.01,
                                        detect_table[,names[k]], results_table[,(index+1):dim(results_table)[2]])
                    colnames(results_table)[((index+1):(index+2))]=c(paste(names[k],".Detected",sep=""),paste(names[k],".Detect_pval",sep=""))
                  }
                  colnames(results_table)[(dim(results_table)[2])]="Significant"
                  index<-which(colnames(results_table)==names[1])
                  if(rpl.mean==TRUE){
                    detection.cols=grep(".Detected",colnames(results_table))
                    Overall.detection=rowSums(results_table[,detection.cols])
                  } else{
                    Overall.detection=detectionCall(subdata.comp, type="probe")
                  }
                  temp<-match(row.names(results_table),names(Overall.detection))
                  Overall.detection<-Overall.detection[temp]
                  results_table=cbind(results_table[,1:(index-1)], Overall.detection, results_table[,index:dim(results_table)[2]])
                  
                  
                  results_table_005<-results_table[results_table$adj.P.Val<0.05,]
                  results_table_001<-results_table[results_table$adj.P.Val<0.01,]
                  results_sum1<-c(paste(name1,name2,sep="_vs_"),dim(exprs(subdata))[1], dim(results_table_001)[1], 
                                  dim(results_table_005)[1])
                  if(is.na(results_sum[1,1])){
                    results_sum[1,]<-results_sum1
                  }else{
                    results_sum<-rbind(results_sum, results_sum1)
                  }
                  
                  
                  volca<-ggplot(data=results_table, aes(x=logFC, y=-log10(adj.P.Val), colour=Significant))+geom_point()+scale_colour_manual(values=c("blue","red"))
                  vplotname=paste("volcano_plot_",levels_groups[i],"_", levels_groups[l],".png", sep="")
                  ggsave(filename = vplotname, path=paste(newDir,"/Graphs", sep=""),scale=1, width=13.5, height=6.68)
                  plotPath=paste("Comparisons/",name1,"_vs_",name2,"/Graphs/",vplotname, sep="")
                  listVolca<-c(listVolca,plotPath)
                  #                   print(listVolca)
                  compName<-c(compName,paste(name1, " vs ", name2,sep=""))
                  
                  resultsFile=paste(newDir,"/",levels_groups[i],"_vs_", levels_groups[l],".csv",sep="")
                  resultsFileName=paste(levels_groups[i],"_vs_", levels_groups[l],".csv",sep="")
                  write.csv(file=resultsFile, col.names=c("ProbeId",colnames(results_table)), results_table)
                  
                  filePath=paste("Comparisons/",name1,"_vs_",name2,"/",resultsFileName,sep="")
                  fileList<-c(fileList,filePath)
                  fileNameList<-c(fileNameList,resultsFileName)
                  
                  # Heatmap
                  if(is.na(heatmap.geneList) | is.null(heatmap.geneList)) {
                    print(colnames(results_table))
                    geneList<-results_table$PROBE_ID[1:100]
                  } else {
                    geneList<-heatmap.geneList
                  }
                  print(geneList)
                  plotname_hm=paste("heatmap_",levels_groups[i],"_", levels_groups[l],".png", sep="")
                  png(filename=paste(newDir,"/Graphs/",plotname_hm, sep=""),  width = 480, height = 780)
                  heatmap_table<-createHeatmap(subdata, geneList=geneList)
                  dev.off()
                  if(arrayType=="Mouse"){
                    heatmap_table<-setAnnotation(as.data.frame(heatmap_table),arrayType = "Mouse")
                  } else {
                    heatmap_table<-setAnnotation(as.data.frame(heatmap_table),arrayType = "Human")
                  }
                  
                  plotPath_hm=paste("Comparisons/",name1,"_vs_",name2,"/Graphs/",plotname_hm, sep="")
                  listHeatMap<-c(listHeatMap,plotPath_hm)
                  
                  resultsFile_hm=paste(newDir,"/",levels_groups[i],"_vs_", levels_groups[l],"_heatmap.csv",sep="")
                  resultsFileName_hm=paste(levels_groups[i],"_vs_", levels_groups[l],"_heatmap.csv",sep="")
                  write.csv(file=resultsFile_hm, heatmap_table)
                  
                  filePath=paste("Comparisons/",name1,"_vs_",name2,"/",resultsFileName_hm,sep="")
                  
                }
              }
              write.table(results_sum, file="comparisons_summary.txt", sep="\t",row.names = FALSE)
            }
          }
)

setGeneric("createHeatmap", function(subdata,geneList=NA) standardGeneric("createHeatmap"))
setMethod("createHeatmap", "IlluminaData", 
          function(subdata, geneList="data.frame"){
            heatmap.data<-exprs(subdata)[geneList,]
            if(require(RColorBrewer)){
              # annotations<-nuID2IlluminaID(featureNames(object), species=organism)
              # positions<-grep("TRUE",geneList[,1] %in% annotations[,5])
              # rownames(exprs(subset))<-annotations[featureNames(subset),4]
              if(require(gplots)){  
                my_palette <- colorRampPalette(c("green","black","red"))(n = 299)
                x.lab.margin=max(nchar(as.character(groups(subdata)[,1])))+(max(nchar(as.character(groups(subdata)[,1])))*0.25)
                out=heatmap.2(heatmap.data, col=my_palette, trace="none", labRow = "", scale="col", margins = c(15, 5), cexCol=1.4)
                carpet_unt<-t(out$carpet)
                row.names(carpet_unt)<-out$rowInd
                heatmap.data_sorted<-heatmap.data[as.numeric(row.names(carpet_unt)[length(row.names(carpet_unt)):1]),]
              }
            }
            return(heatmap.data_sorted)
          }
)

setGeneric("removeSamples", function(object,remove) standardGeneric("removeSamples"))
setMethod("removeSamples", "IlluminaData",
          function(object, remove=""){
            sample.names=sampleNames(object)
            sample.index=which(sample.names %in% remove)
            object=object[,-sample.index]
            object@groups<-groups(object)[-sample.index,]
            object@groups<-droplevels(object@groups)
            return(object)
          }
)


setGeneric("compareGroup", function(object, pairedTest= FALSE) standardGeneric("compareGroup"))
setMethod("compareGroup", "IlluminaData",
          function(object, pairedTest){
            if(all(is.na(object@designMatrix))){
              print("IlluminaData object has no designMatrix defined, use createDesign to create one") 
            
            }else{
              
              design<-object@designMatrix
              fit <- lmFit(object, design)
              first_cont=1
              for(i in 1:dim(design)[2]){
                if(i!=dim(design)[2]){
                  z=i+1
                  for(l in z:dim(design)[2]){
                    if(first_cont==1){
                      first_cont=0
                      title=paste(as.character(colnames(design)[i]),"vs",as.character(colnames(design)[l]),sep="")
                      contrast=paste(as.character(colnames(design)[i]),"-",as.character(colnames(design)[l]),sep="") 
                      global_contrasts=contrast
                    }else{
                      title=paste(as.character(colnames(design)[i]),"vs",as.character(colnames(design)[l]),sep="")
                      contrast=paste(as.character(colnames(design)[i]),"-",as.character(colnames(design)[l]),sep="") 
                      global_contrasts=c(global_contrasts,contrast)
                    }
                  }
                  contrast.matrix<-makeContrasts(contrasts=global_contrasts, levels=colnames(design))
                }
              }  
              
              fit2<-contrasts.fit(fit,contrast.matrix)
              fit2<-eBayes(fit2)
              ncontrasts=dim(contrast.matrix)[2]
              if (pairedTest==TRUE){
                ncontrasts=1
                print("paired")
              }
              
              first=1
              for(i in 1:ncontrasts){
                if(first==1){
                  first=0
                  resultsData<-topTable(fit2, adjust="BH", number=dim(object)[1], coef=i)
                  colnames(resultsData)[c(10,14)]<-paste(colnames(resultsData[,c(10,14)]),colnames(contrast.matrix)[i],sep=".")
                #  results<-resultsData[order(resultsData[,1])]
                }else{
                  top<-topTable(fit2, adjust="BH", number=dim(object)[1],coef=i)
                  colnames(top)<-paste(colnames(top),colnames(contrast.matrix)[i], sep=".")
                  top<-top[order(top[,1]),]
                  resultsData<-cbind(resultsData,top[,c(11,14)])
                }
              }
            }
            exprs_data<-exprs(object)
            sorted_exprs<-exprs_data[match(resultsData[,1],row.names(exprs_data)),]
            resultsData<-cbind(resultsData,sorted_exprs)
           # resultsData$Significant<-as.factor(abs(resultsData$logFC)>1.5 & resultsData$adj.P.Val<0.05)
            resultsData$Significant<-resultsData$adj.P.Val<0.05
            results(object)<-resultsData
            return(object)
            
  }
)

setGeneric("createDesign", function(object,group=NA,pairedTest=FALSE) standardGeneric("createDesign"))
setMethod("createDesign", "IlluminaData",
          function(object,group, pairedTest){
            if(is.na(group)){
              group<-groups(object)
            } 
            if (pairedTest==TRUE){
              SibShip <- factor(group$Sibship)
              Group <- factor(group$Group)
              designMatrix <- model.matrix(~0+Group+SibShip)
              colnames(designMatrix)[1:2]<-levels(as.factor(group$Group))
              designMatrix(object)<-designMatrix
              print("paired")
            } else {
              designMatrix<-model.matrix(~0+as.character(group$Group))
              colnames(designMatrix)<-levels(as.factor(group$Group))
              designMatrix(object)<-designMatrix
            }
            return(object)
})

setGeneric("setControlData", function(object, filename) standardGeneric("setControlData"))
setMethod("setControlData","IlluminaData",function(object,filename="character"){
  controlData<-getControlData(filename)
  controlData(object)<-controlData
  return(object)
})


### No replicate analysis
setGeneric("multipleComparisonsNoRpl", function(object, arrayType, heatmap.geneList=NA) standardGeneric("multipleComparisonsNoRpl"))
setMethod("multipleComparisonsNoRpl", "IlluminaData",
          function(object, arrayType="", heatmap.geneList=NA){
            if(is.na(groups(object))){
              print("Groups is empty, needed to do any comparisons")
            }else{
              #png()
              results_sum<-c("Comparisons","Number of selected probes","logFC>=1", "logFC>=2")
              listVolca<-vector(mode="list")
              listHeatMap<-vector(mode="list")
              compName<-vector(mode="list")
              fileNameList<-vector(mode="list")
              fileList<-vector(mode="list")
              fileNameList_hm<-vector(mode="list")
              fileList_hm<-vector(mode="list")
              levels_groups=levels(as.factor(groups(object)[,2]))
              group_all=groups(object)
              for(i in 1:(length(levels_groups)-1)){
                for(l in (i+1):length(levels_groups)){
                  #create folder
                  name1=levels_groups[i]
                  name2=levels_groups[l]
                  newDir=file.path(getwd(),paste("Comparisons/",name1,"_vs_",name2,sep=""))
                  #print(newDir)
                  dir.create(newDir,recursive=T)
                  dir.create(paste(newDir,"/Graphs", sep=""))
                  
                  subgroup<-group_all[group_all[,2]==levels_groups[i] | group_all[,2]==levels_groups[l],]
                  subdata<-processData(object,subgroup)
                  subdata.comp=cbind(exprs(subdata)[,2]-exprs(subdata[,1]),exprs(subdata))  # neg=downregulates in second sample
                  colnames(subdata.comp)[1]="logFC"
                  results_table<-as.data.frame(subdata.comp)
                  results_table<-results_table[order(results_table$logFC),]
                  results_table<-setAnnotation(as.data.frame(results_table), arrayType = "Mouse")
                  results_table=cbind(PROBE_ID=rownames(results_table),results_table)
                  
                  temp=match(row.names(results_table),row.names(detection(subdata)))
                  detect_table=as.data.frame(detection(subdata)[temp,])
                  names=colnames(detect_table)
                  for (k in 1:length(names)){ 
                    index<-which(colnames(results_table)==names[k])
                    results_table=cbind(results_table[,1:index], detect_table[,names[k]]<=0.01,
                                        detect_table[,names[k]], results_table[,(index+1):dim(results_table)[2]])
                    colnames(results_table)[((index+1):(index+2))]=c(paste(names[k],".Detected",sep=""),paste(names[k],".Detect_pval",sep=""))
                  }
                  colnames(results_table)[(dim(results_table)[2])]="Symbol"
                  index<-which(colnames(results_table)==names[1])
                  Overall.detection=detectionCall(subdata, type="probe") # subdata or subdata.comp?
                  temp<-match(row.names(results_table),names(Overall.detection))
                  Overall.detection<-Overall.detection[temp]
                  results_table=cbind(results_table[,1:(index-1)], Overall.detection, results_table[,index:dim(results_table)[2]])
                  
                  results_table_logFC1<-results_table[abs(results_table$logFC)>=1,]
                  results_table_logFC3<-results_table[abs(results_table$logFC)>=2,]
                  results_sum1<-c(paste(name1,name2,sep="_vs_"),dim(subdata.comp)[1], dim(results_table_logFC1)[1], 
                                  dim(results_table_logFC3)[1])
                  results_sum<-rbind(results_sum, results_sum1)
                  
                  plotname=paste("logFC_hist",levels_groups[i],"_", levels_groups[l],".png", sep="")
                  png(filename=paste(newDir,"/Graphs/",plotname, sep=""))
                  hist<-hist(results_table$logFC, breaks=10000, width = 880, height = 880, xlab="logFC", main="")
                  dev.off()
                  plotPath=paste("Comparisons/",name1,"_vs_",name2,"/Graphs/",plotname, sep="")
                  listVolca<-c(listVolca,plotPath)
                  #                   print(listVolca)
                  compName<-c(compName,paste(name1, " vs ", name2,sep=""))
                  
                  resultsFile=paste(newDir,"/",levels_groups[i],"_vs_", levels_groups[l],".csv",sep="")
                  resultsFileName=paste(levels_groups[i],"_vs_", levels_groups[l],".csv",sep="")
                  write.csv(file=resultsFile, col.names=c("ProbeId",colnames(results_table)), results_table)
                  #write.csv(file=resultsFile, results_table)
                  filePath=paste("Comparisons/",name1,"_vs_",name2,"/",resultsFileName,sep="")
                  fileList<-c(fileList,filePath)
                  fileNameList<-c(fileNameList,resultsFileName)
                  
                  # Heatmap
                  if(is.na(heatmap.geneList)) {
                    geneList<-results_table$PROBE_ID[1:100]
                  } else {
                    geneList<-heatmap.geneList
                  }
                  plotname_hm=paste("heatmap_",levels_groups[i],"_", levels_groups[l],".png", sep="")
                  png(filename=paste(newDir,"/Graphs/",plotname_hm, sep=""))
                  heatmap_table<-createHeatmap(subdata, geneList=geneList)
                  dev.off()
                  if(arrayType=="Mouse"){
                    heatmap_table<-setAnnotation(as.data.frame(heatmap_table),arrayType = "Mouse")
                  } else {
                    heatmap_table<-setAnnotation(as.data.frame(heatmap_table),arrayType = "Human")
                  }
                  
                  plotPath_hm=paste("Comparisons/",name1,"_vs_",name2,"/Graphs/",plotname_hm, sep="")
                  listHeatMap<-c(listHeatMap,plotPath_hm)
                  
                  resultsFile_hm=paste(newDir,"/",levels_groups[i],"_vs_", levels_groups[l],"_heatmap.csv",sep="")
                  resultsFileName_hm=paste(levels_groups[i],"_vs_", levels_groups[l],"_heatmap.csv",sep="")
                  write.csv(file=resultsFile_hm, heatmap_table)
                  
                  filePath=paste("Comparisons/",name1,"_vs_",name2,"/",resultsFileName_hm,sep="")
                  fileList_hm<-c(fileList_hm,filePath)
                  fileNameList_hm<-c(fileNameList_hm,resultsFileName_hm)
                }
              }
              reportName="Multiple Comparisons results summary"
              htmlResults<-HTMLReport(shortName=reportName, title="Results Report",center=T, reportDirectory=".")
              publish(hwrite(reportName, heading=3, center=T),htmlResults)
              publish(hwrite(results_sum,table=T,row.names=F,col.names=T, center=T), htmlResults)
              #print(length(listVolca[]))
              # print(compName[[1]][1])
              for(i in 1:(length(listVolca[]))){
                # print(i)
                publish(hwrite(as.character(compName[[i]][1]), heading=3, center=T),htmlResults) 
                #print(listVolca[[i]][1])
                publish(hwrite("LogFC distribution", center=F, style="font-weight:bold;  font-size:120%"), htmlResults)
                img<-hwriteImage(listVolca[[i]][1], link=listVolca[[i]][1],image.border=4,center=T, width=100*6, height=100*6)
                publish(hwrite(img,br=T), htmlResults)
                #textVolcano="Volcano Plot, significance as a function of the log fold change, the adjusted p-value is transformed to a negative log10 scale
                #and plotted on the y axis, results are significant if above 1.3 (adj.P.Val<0.05), 2 for higher significance (adj.P.val<0.01)."
                #publish(hwrite(textVolcano, br=F, center=F), htmlResults)
                textCompa1="Link to access the results table:"
                #textCompa2="To sort your data use the logFC column (log2 scale) and adjusted pvalue (adj.P.val), the p-value should be at least below 0.05,
                #the fold change can be set according to the magnitude of the change you are looking for."
                publish(hwrite(textCompa1, br=F, center=F), htmlResults)
                publish(hwrite(fileNameList[[i]][1],link=fileList[[i]][1],center=F,br=T),htmlResults)
                #publish(hwrite(textCompa2, br=F, center=F), htmlResults) 
                
                publish(hwrite("space", br=F, center=F, style="color:white"), htmlResults) 
                publish(hwrite("HEATMAP", center=F, style="font-weight:bold;  font-size:120%"), htmlResults)
                publish(hwrite(" ", heading=3, center=T),htmlResults) 
                img<-hwriteImage(listHeatMap[[i]][1], link=listHeatMap[[i]][1], width=580, image.border=4, center=T)
                publish(hwrite(img,br=T), htmlResults)
                #textHeatMap="Heatmap plot. Top 100 genes sorted by adj.P.Val were grouped using hierarchical clustering. Note that not all of the 100 genes are necessarily significant. 
                #Intensity values were centered and scaled across samples. Green color indicates lower intensity and Red color indicates higher intensity. A light blue traceback line in the color key indicates the distribution 
                #of scaled intensity values."
                #publish(hwrite(textHeatMap, br=F, center=F), htmlResults)
                publish(hwrite("Link to access the heatmap table", br=F, center=F), htmlResults)
                publish(hwrite(fileNameList_hm[[i]][1],link=fileList_hm[[i]][1],center=F,br=T),htmlResults)
                textHeatMap_table="This table contains intensity values and genes are in the same order as presented in the heatmap plot."
                publish(hwrite(textHeatMap_table, br=F, center=F), htmlResults)
                
              }
              finish(htmlResults)
              #publish(hwrite(paste("Group ",group_cat[i]," ",group_cat[l], " Scatter plot", sep=""), heading=2, center=T),htmlControl)
              #img<-hwriteImage(paste(folderGraph2,"/",nameOnly,sep=""), link=paste(folderGraph2,"/",nameOnly,sep=""),center=T)
              #publish(hwrite(img,br=T), htmlControl)
            }
          }
)
#####

setGeneric("results",function(object) standardGeneric("results"))
setMethod("results", "IlluminaData", function(object) object@results)

setGeneric("groups",function(object) standardGeneric("groups"))
setMethod("groups", "IlluminaData", function(object) object@groups)

setGeneric("designMatrix",function(object) standardGeneric("designMatrix"))
setMethod("designMatrix", "IlluminaData", function(object) object@designMatrix)

setGeneric("controlData",function(object) standardGeneric("controlData"))
setMethod("controlData", "IlluminaData", function(object) object@controlData)


setGeneric("designMatrix<-", function(object, value) standardGeneric("designMatrix<-"))
setReplaceMethod("designMatrix", "IlluminaData", function(object, value){
  object@designMatrix=value
  object
})

setGeneric("results<-", function(object, value) standardGeneric("results<-"))
setReplaceMethod("results", "IlluminaData", function(object, value){
  object@results=value
  object
})

setGeneric("controlData<-", function(object, value) standardGeneric("controlData<-"))
setReplaceMethod("controlData", "IlluminaData", function(object, value){
  object@controlData=value
  object
})


######LOG removed:
#v1.01 : removed the QC report and biological report that generated html results
#v1.01: remove the mutliple comparisons method generated html report