rm(list=ls())
require(dplyr)
require(tidyr)
require(zetadiv)
require(geosphere)
require(relaimpo)
require(tidyverse)

wd <- "/staging/sn1/alsimons/CoastalBMIs" #On the cluster
wd <- "~/Desktop/CoastalBMIs/" #Local drive
setwd(wd)

#Read in taxa abundance by sample.
communityInput <- read.table("Master - Benthos ed12.csv", header=T, sep=",",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")

#Create unique ID per sample merging station ID with sample replicate.
communityInput <- unite(communityInput,SampleID,c("StationID","Replicate"),sep="_",remove=FALSE)

#Create basic sample metadata names list.
sampleNames <- c("SampleID","StationWaterDepth","latitude","longitude")

#Remove empty taxon rows.
communityInput <- communityInput[communityInput$Taxon!="",]

#Remove duplicate entries.
communityInput <- communityInput[!duplicated(communityInput[,c("SampleID","Taxon")]),]

#Read in sample station metadata.
stationInput <- read.table("Master Data - Station Info.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8",allowEscapes=TRUE)

#Merge in sample sation metadata to community data.
communityInput <- dplyr::left_join(communityInput,stationInput,by=c("StationID"))

#Read in BMI taxonomic data.
taxaInput <- read.table("Master - Taxonomic Info edit.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8",quote="")
#Remove quote characters.
taxaInput <- as.data.frame(sapply(taxaInput, function(x) gsub("\"", "", x)))

#Merge in taxonomic data into community data.
communityInput <- dplyr::left_join(communityInput,taxaInput,by=c("Taxon"))

#Read in assessment scores.
assessmentInput <- read.table("Master Data - Assessment Scores.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Merge in assessment scores into community data.
communityInput <- dplyr::left_join(communityInput,assessmentInput,by=c("StationID","Replicate"))

#Read in sediment data.
sedimentInput <- read.table("Master Data - Grainsize.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Reformat sediment data so the percent of each sediment type is a column.
sedimentInput <- sedimentInput %>% group_by(GrainsizeBin) %>% mutate(n = 1:n()) %>% tidyr::spread(GrainsizeBin,Percent)
sedimentInput$n <- NULL
#Get sediment type names.
sedimentNames <- colnames(sedimentInput)
sedimentNames <- sedimentNames[sedimentNames != "StationID"]

#Merge in sediment data into community data.
communityInput <- dplyr::left_join(communityInput,sedimentInput,by=c("StationID"))

#Read in chemical composition data.
chemInput <- read.table("Data - chem un-edited.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Get the percent occurrence of each chemical analyte.
chemFreq <- as.data.frame(table(chemInput$Parameter)/length(unique(chemInput$StationID)))
colnames(chemFreq) <- c("Parameter","Freq")
chemFreq$Parameter <- as.character(chemFreq$Parameter)

#Select the most commonly occurring chemical analytes.
chemCutoff <- 0.70
chemList <- chemFreq[chemFreq$Freq >= chemCutoff, "Parameter"]
chemSubset <- chemInput[chemInput$Parameter %in% chemList,]

#Reformat chemical data so the concentration of each analyte is a column.
chemSubset <- chemSubset[,c("StationID","Parameter","Result")]
chemSubset <- chemSubset %>% group_by(Parameter) %>% mutate(n = 1:n()) %>% tidyr::spread(Parameter,Result)
chemSubset$n <- NULL

#Convert negative concentrations to NA
chemSubset[chemSubset < 0] <- NA

#Calculate the mean analyte concentration per sample.
chemSubset <- aggregate(. ~ StationID, chemSubset, function(x) mean(x, na.rm=T), na.action = na.pass)
#Convert NaN to NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
chemSubset[is.nan(chemSubset)] <- NA

#Get chemical names.
chemNames <- colnames(chemSubset)
chemNames <- chemNames[chemNames != "StationID"]

#Merge in chemical composition data into community data.
communityInput <- dplyr::left_join(communityInput,chemSubset,by=c("StationID"))

#Choose a taxonomic level to group count data by.
#Levels are Domain, Kingdom, Phylum, Class, Order, Family, GenusSpecies, OTUID
taxonomicLevels <- c("Phylum","Class","Order","Family","Genus","species")
taxonomicLevel <- c("Order") #Choose a taxonomic level to aggregate count data on.

#Create simplified community data frame to summarize data at a particular taxonomic level.
communityInputSummarized <- communityInput[, c("SampleID",taxonomicLevel,"Abundance")]
#Remove rows without assigned taxonomy at a particular level.
communityInputSummarized <- communityInputSummarized[communityInputSummarized[,c(taxonomicLevel)]!="",]
#Aggregate communities at a particular taxonomic level.
communityInputSummarized <- as.data.frame(aggregate(formula(paste0("Abundance ~ SampleID+",taxonomicLevel)),communityInputSummarized,sum,na.action = na.omit))

#Create master taxa list.
uniqueTaxa <- as.data.frame(unique(communityInputSummarized[,c(taxonomicLevel)]))
colnames(uniqueTaxa) <- c("Taxon") #Rename column
uniqueTaxa$Taxon <- as.character(uniqueTaxa$Taxon) #Coerce data to character type from factor.
uniqueTaxa <- dplyr::rename(uniqueTaxa, !!taxonomicLevel := Taxon) #Rename taxonomic level to just 'Taxon'.

#Create a community matrix where rows are the taxa IDs and the columns are unique sampleIDs
communityPA <- uniqueTaxa
for(sample in unique(communityInputSummarized$SampleID)){
  tmp <- communityInputSummarized[communityInputSummarized$SampleID==sample,c("SampleID",taxonomicLevel,"Abundance")] #Subset community data by sample.
  tmp <- merge(tmp,uniqueTaxa,by=taxonomicLevel,all=TRUE)
  rownames(tmp) <- tmp$Taxon #Assign taxa ID to row names
  tmp <- dplyr::rename(tmp, !!sample := Abundance) #Assign sampled ID to the taxa abundance present in that sample.
  tmp <- as.data.frame(tmp[, -which(names(tmp) %in% c("SampleID",taxonomicLevel))])
  colnames(tmp) <- sample
  communityPA <- cbind(communityPA,tmp)
}
#Convert abundance data to presence/absence.
rownames(communityPA) <- communityPA[,c(taxonomicLevel)] #Assign taxa ID to row names
communityPA <- as.data.frame(communityPA[, -which(names(communityPA) %in% c(taxonomicLevel))])
#Convert community matrix to presence/absence
communityPA[is.na(communityPA)] <- 0
communityPA[communityPA >= 1] <- 1

#Group samples by stratum / year / condition score.  Then keep groups with at least a certain number of samples.
sampleList <- communityInput[,c("SampleID","Stratum","SampleYear","ConditionScore")]
sampleList <- sampleList[!duplicated(sampleList),]
sampleList <- sampleList[!is.na(sampleList$ConditionScore),]
tmp <- as.data.frame(table(sampleList[,c("Stratum","SampleYear","ConditionScore")]))
tmp$SampleYear <- as.numeric(as.character(tmp$SampleYear))
tmp$ConditionScore <- as.numeric(as.character(tmp$ConditionScore))
tmp$Stratum <- as.character(tmp$Stratum)
sampleList <- dplyr::left_join(sampleList,tmp,by=c("Stratum","SampleYear","ConditionScore"))
sampleMax <- 15 #Minimum number of samples per stratum / year / condition score group to be considered for analysis.
sampleList <- sampleList[sampleList$Freq >= sampleMax,]

zetaMax <- 10 #Maximum order to calculate zeta diversity.
sampleNum <- 12 #Number of samples used at a time to perform zeta diversity analysis on.
zetaAnalysis <- data.frame()
for(j in 1:100){
  for(stratumType in unique(sampleList$Stratum)){
    for(year in unique(sampleList$SampleYear)){
      for(score in unique(sampleList$ConditionScore)){
        subList <- sampleList[sampleList$Stratum==stratumType & sampleList$SampleYear==year & sampleList$ConditionScore==score,]
        if(nrow(subList) > 0){
          sampleCommunity <- communityPA[,colnames(communityPA) %in% subList$SampleID]
          data.spec <- as.data.frame(t(sampleCommunity))
          data.spec <- data.spec[sample(nrow(data.spec),sampleNum),]
          zetaDecay <- Zeta.decline.ex(data.spec,orders=1:zetaMax,rescale=TRUE,plot=FALSE)
          zeta_1 <- Zeta.order.ex(data.spec,order=1,rescale=TRUE)$zeta.val #lower order zeta diversity measure.
          zeta_1sd <- Zeta.order.ex(data.spec,order=1,rescale=TRUE)$zeta.val.sd #lower order zeta diversity measure standard deviation.
          zeta_2 <- zetaDecay$zeta.val[2] #lower order zeta diversity measure.
          zeta_2sd <- zetaDecay$zeta.val.sd[2] #lower order zeta diversity measure standard deviation.
          zeta_3 <- zetaDecay$zeta.val[3] #lower order zeta diversity measure.
          zeta_3sd <- zetaDecay$zeta.val.sd[3] #lower order zeta diversity measure standard deviation.
          zeta_4 <- zetaDecay$zeta.val[4] #lower order zeta diversity measure.
          zeta_4sd <- zetaDecay$zeta.val.sd[4] #lower order zeta diversity measure standard deviation.
          zeta_N <- zetaDecay$zeta.val[zetaMax] #Higher order zeta diversity measure.
          zeta_Nsd <- zetaDecay$zeta.val.sd[zetaMax] #Higher order zeta diversity measure standard deviation.
          ExpExp <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
          ExpAIC <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
          PLExp <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
          PLAIC <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
          metadataSubset <- communityInput[communityInput$SampleID %in% rownames(data.spec),]
          metadataSubset <- metadataSubset[,colnames(metadataSubset) %in% c(sampleNames,chemNames,sedimentNames)]
          metadataSubset <- metadataSubset[!duplicated(metadataSubset),]
          meanDist <- mean(distm(metadataSubset[,c("longitude","latitude")]))
          sdDist <- sd(distm(metadataSubset[,c("longitude","latitude")]))
          meanDepth <- mean(metadataSubset$StationWaterDepth)
          sdDepth <- sd(metadataSubset$StationWaterDepth)
          chemData <- metadataSubset[,colnames(metadataSubset) %in% chemNames]
          chemRow <- colMeans(chemData,na.rm=T)
          chemRow[is.nan(chemRow)] <- NA
          sedimentData <- metadataSubset[,colnames(metadataSubset) %in% sedimentNames]
          sedimentRow <- colMeans(sedimentData,na.rm=T)
          sedimentRow[is.nan(sedimentRow)] <- NA
          dataRow <- t(cbind(t(as.data.frame(list(c(stratumType,year,score,meanDepth,sdDepth,meanDist,sdDist,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,zeta_2,zeta_2sd,zeta_3,zeta_3sd,zeta_4,zeta_4sd,ExpExp,ExpAIC,PLExp,PLAIC,as.list(as.numeric(sedimentRow)),as.list(as.numeric(chemRow))))))))
          colnames(dataRow) <- c("stratumType","year","ConditionScore","meanDepth","sdDepth","meanDist","sdDist","zeta_N","zeta_Nsd","zeta_1","zeta_1sd","zeta_2","zeta_2sd","zeta_3","zeta_3sd","zeta_4","zeta_4sd","ExpExp","ExpAIC","PLExp","PLAIC",sedimentNames,chemNames)
          rownames(dataRow) <- NULL
          zetaAnalysis <- rbind(zetaAnalysis,dataRow)
          print(paste(j,stratumType,year,score,meanDepth,sdDepth,meanDist,sdDist,zeta_N,zeta_Nsd,zeta_1,zeta_1sd,zeta_2,zeta_2sd,zeta_3,zeta_3sd,zeta_4,zeta_4sd,ExpExp,ExpAIC,PLExp,PLAIC))
          print(chemRow)
          print(sedimentRow)
        }
      }
    }
  }
}

colnames(zetaAnalysis) <- c("stratumType","year","ConditionScore","meanDepth","sdDepth","meanDist","sdDist","zeta_N","zeta_Nsd","zeta_1","zeta_1sd","zeta_2","zeta_2sd","zeta_3","zeta_3sd","zeta_4","zeta_4sd","ExpExp","ExpAIC","PLExp","PLAIC",sedimentNames,chemNames)
tmp <- zetaAnalysis$stratumType
indx <- sapply(zetaAnalysis, is.factor)
zetaAnalysis[indx] <- lapply(zetaAnalysis[indx], function(x) as.numeric(as.character(x)))
zetaAnalysis <- do.call(data.frame,lapply(zetaAnalysis, function(x) replace(x, is.infinite(x),NA)))
zetaAnalysis$stratumType <- tmp
#Save zeta diversity analysis for a given taxonomic level.
write.table(zetaAnalysis,paste("coastalBMIs",taxonomicLevel,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#

##To run locally.
zetaAnalysis <- read.table(paste("coastalBMIs",taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

zetaModel <- lm(ConditionScore~zeta_1+zeta_2+zeta_N,data=zetaAnalysis)
plot(zetaModel$model$ConditionScore,zetaModel$fitted.values)
cor.test(zetaModel$model$ConditionScore,zetaModel$fitted.values)
zetaAnalysis$modeledConditionScore <- zetaModel$fitted.values
calc.relimp(zetaModel)
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(zetaModel)
dev.off()

# Assessing R2 shrinkage using 10-Fold Cross-Validation 
require(bootstrap)
require(caret)
set.seed(1)
train.control <- trainControl(method="repeatedcv",number=10,repeats=10)
ConditionScoremodel <- train(ConditionScore~zeta_1+zeta_2+zeta_N,data=zetaAnalysis,method="lm",trControl=train.control)
print(ConditionScoremodel)

require(ggplot2)
require(viridis)
zetaPlot <- ggplot(zetaAnalysis, aes(x=ConditionScore,y=modeledConditionScore,color=meanDepth))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=modeledConditionScore))
zetaPlot+xlab("Condition Score")+ylab("Modeled Condition Score")+scale_color_gradientn("Mean Depth (m)",colours = rev(plasma(10)))

#Check for correlation patterns between zeta diversity and environmental parameters.
require("PerformanceAnalytics")
chart.Correlation(zetaAnalysis[,c("modeledConditionScore","ConditionScore","meanDepth","meanDist","zeta_1","zeta_2","zeta_3","zeta_4","zeta_N")], histogram=TRUE, method="pearson")

#Mapping metadata
require(ggplot2)
require(ggmap)
require(maps)
require(mapview)
require(mapdata)
require(munsell)
require(leaflet)
require(devtools)
require(webshot)
require(viridis)
#Map data for continuous variables.
CalMap = leaflet(stationInput) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=viridis(10),domain=stationInput$StationWaterDepth)
CalMap %>% addCircleMarkers(color = ~ColorScale(StationWaterDepth), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~StationWaterDepth,title="SCCWRP sample<br>Depth (m)")

#Map data for categorical variables.
pal <- colorFactor(palette = 'viridis', domain = stationInput$Stratum)
CalMap = leaflet(stationInput) %>% 
  addTiles()
CalMap %>% addCircleMarkers(color = ~pal(Stratum), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=pal,values=~Stratum,title="SCCWRP sample<br>Regions")
