rm(list=ls())
require(dplyr)
require(tidyr)
require(zetadiv)
require(geosphere)
require(relaimpo)
require(tidyverse)

wd <- "/scratch/alsimons/CoastalBMIs" #On the cluster
wd <- "~/Desktop/CoastalBMIs/" #Local drive
setwd(wd)

#Read in taxa abundance by sample.
communityInput <- read.table("Master - Benthos ed12.csv", header=T, sep=",",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")

#Create unique ID per sample merging station ID with sample replicate.
communityInput <- unite(communityInput,SampleID,c("StationID","Replicate"),sep="_",remove=FALSE)

#Create basic sample metadata names list.
sampleNames <- c("SampleID","StationWaterDepth","DistanceToOutfall","latitude","longitude")

#Remove empty taxon rows.
communityInput <- communityInput[communityInput$Taxon!="",]

#Remove duplicate entries.
communityInput <- communityInput[!duplicated(communityInput[,c("SampleID","Taxon")]),]

#Read in sample station metadata.
stationInput <- read.table("Master Data - Station Info.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8",allowEscapes=TRUE)

#Merge in sample sation metadata to community data.
communityInput <- dplyr::left_join(communityInput,stationInput,by=c("StationID"))

#Read in distance to outfall pipe per sample location.
outfallInput <- read.table("SampleStationsWithDistanceToOutall.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8",allowEscapes=TRUE)
#Keep only relevant columns.
outfallInput <- outfallInput[,c("StationID","DistanceToOutfall")]

#Merge in distance to outfall pipe locations into the community data.
communityInput <- dplyr::left_join(communityInput,outfallInput,by=c("StationID"))

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

#Get the percent occurrence o each sediment type.
sedimentFreq <- as.data.frame(table(sedimentInput$GrainsizeBin)/length(unique(sedimentInput$StationID)))
colnames(sedimentFreq) <- c("GrainsizeBin","Freq")
sedimentFreq$GrainsizeBin <- as.character(sedimentFreq$GrainsizeBin)

#Select the most commonly occurring sediment types.
sedimentCutoff <- 0.70
sedimentList <- sedimentFreq[sedimentFreq$Freq >= sedimentCutoff, "GrainsizeBin"]
sedimentInput <- sedimentInput[sedimentInput$GrainsizeBin %in% sedimentList,]

#Reformat sediment data so the percent of each sediment type is a column.
sedimentInput <- sedimentInput %>% dplyr::group_by(GrainsizeBin) %>% dplyr::mutate(n = 1:n()) %>% tidyr::spread(GrainsizeBin,Percent)
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
chemSubset <- chemSubset %>% dplyr::group_by(Parameter) %>% dplyr::mutate(n = 1:n()) %>% tidyr::spread(Parameter,Result)
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

set.seed(1)
zetaMax <- 10 #Maximum order to calculate zeta diversity.
zetaNames <- paste("zeta",seq(1:zetaMax),sep="_") #Generate list of names for zeta diversity measures.
zetaSDNames <- paste("zeta",seq(1:zetaMax),"sd",sep="_") #Generate list of names for the standard deviation on zeta diversity measures.
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
          ExpExp <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
          ExpAIC <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
          PLExp <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
          PLAIC <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
          zetaOrders <- zetaDecay$zeta.val #Generate list of zeta diversity measures from 1 to zetaMax.
          zetaOrders[1] <- Zeta.order.ex(data.spec,order=1,rescale=TRUE)$zeta.val #Recalculate zeta_1 as the average taxonomic richness.
          zetaSDOrders <- zetaDecay$zeta.val.sd #Generate list of the standard deviation on zeta diversity measures from 1 to zetaMax.
          zetaSDOrders[1] <- Zeta.order.ex(data.spec,order=1,rescale=TRUE)$zeta.val.sd #Recalculate zeta_1sd as the standard deviation on the average taxonomic richness.
          metadataSubset <- communityInput[communityInput$SampleID %in% rownames(data.spec),]
          metadataSubset <- metadataSubset[,colnames(metadataSubset) %in% c(sampleNames,chemNames,sedimentNames)]
          metadataSubset <- metadataSubset[!duplicated(metadataSubset),]
          meanDist <- mean(distm(metadataSubset[,c("longitude","latitude")]))
          sdDist <- sd(distm(metadataSubset[,c("longitude","latitude")]))
          meanOutfallDist <- mean(metadataSubset$DistanceToOutfall)
          sdOutfallDist <- sd(metadataSubset$DistanceToOutfall)
          meanDepth <- mean(metadataSubset$StationWaterDepth)
          sdDepth <- sd(metadataSubset$StationWaterDepth)
          chemData <- metadataSubset[,colnames(metadataSubset) %in% chemNames]
          chemRow <- colMeans(chemData,na.rm=T)
          chemRow[is.nan(chemRow)] <- NA
          sedimentData <- metadataSubset[,colnames(metadataSubset) %in% sedimentNames]
          sedimentRow <- colMeans(sedimentData,na.rm=T)
          sedimentRow[is.nan(sedimentRow)] <- NA
          dataRow <- t(cbind(t(as.data.frame(list(c(stratumType,year,score,meanDepth,sdDepth,meanDist,sdDist,meanOutfallDist,sdOutfallDist,zetaOrders,zetaSDOrders,ExpExp,ExpAIC,PLExp,PLAIC,as.list(as.numeric(sedimentRow)),as.list(as.numeric(chemRow))))))))
          colnames(dataRow) <- c("stratumType","year","ConditionScore","meanDepth","sdDepth","meanDist","sdDist","meanOutfallDist","sdOutfallDist",zetaNames,zetaSDNames,"ExpExp","ExpAIC","PLExp","PLAIC",sedimentNames,chemNames)
          rownames(dataRow) <- NULL
          zetaAnalysis <- rbind(zetaAnalysis,dataRow)
          print(paste(j,stratumType,year,score,meanDepth,sdDepth,meanDist,sdDist,meanOutfallDist,sdOutfallDist,zetaOrders,zetaSDOrders,ExpExp,ExpAIC,PLExp,PLAIC))
          print(chemRow)
          print(sedimentRow)
        }
      }
    }
  }
}

colnames(zetaAnalysis) <- c("stratumType","year","ConditionScore","meanDepth","sdDepth","meanDist","sdDist","meanOutfallDist","sdOutfallDist",zetaNames,zetaSDNames,"ExpExp","ExpAIC","PLExp","PLAIC",sedimentNames,chemNames)
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

#Model the condition score as a function of zeta diversity measures.
zetaModel <- lm(ConditionScore ~ zeta_1+zeta_2+zeta_10,data=zetaAnalysis)
calc.relimp(zetaModel)
plot(zetaModel$model$ConditionScore,zetaModel$fitted.values)
cor.test(zetaModel$model$ConditionScore,zetaModel$fitted.values)
zetaAnalysis$modeledConditionScore <- zetaModel$fitted.values
#Plot model errors.
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(zetaModel)
dev.off()

#Model the modeled condition score as a function of environmental metrics.
require(car)
environmentNames <- c(chemNames,sedimentNames,"meanDist","meanDepth","meanOutfallDist")
zetaModel <- lm(formula(paste("modeledConditionScore ~ ",paste(environmentNames,collapse="+"))),data=zetaAnalysis)
calc.relimp(zetaModel)
#Filter our variables with a high variance inflation factor.
varList <- car::vif(zetaModel)
varList <- names(varList[varList <= 10])
#Model the modeled condition score as a function of the remaining environmental metrics
zetaModel <- lm(formula(paste("modeledConditionScore ~ ",paste(varList,collapse="+"))),data=zetaAnalysis)
calc.relimp(zetaModel)
cor.test(zetaModel$model$modeledConditionScore,zetaModel$fitted.values)

#Model the condition score as a function of environmental metrics.
zetaModel <- lm(formula(paste("ConditionScore ~ ",paste(environmentNames,collapse="+"))),data=zetaAnalysis)
calc.relimp(zetaModel)
#Filter our variables with a high variance inflation factor.
varList <- car::vif(zetaModel)
varList <- names(varList[varList <= 10])
#Model the condition score as a function of the remaining environmental metrics
zetaModel <- lm(formula(paste("ConditionScore ~ ",paste(varList,collapse="+"))),data=zetaAnalysis)
calc.relimp(zetaModel)
cor.test(zetaModel$model$ConditionScore,zetaModel$fitted.values)

# Assessing R2 shrinkage using 10-Fold Cross-Validation 
require(bootstrap)
require(caret)
set.seed(1)
train.control <- trainControl(method="repeatedcv",number=10,repeats=10)
ConditionScoremodel <- train(ConditionScore~zeta_1+zeta_2+zeta_10,data=zetaAnalysis,method="lm",trControl=train.control)
print(ConditionScoremodel)

#Correlation plots of mean and modeled condition scores against environmental parameters
#with a variance inflation factor under 10 and various measures of zeta diversity.
require(Hmisc)
require(corrplot)
taxonomicLevel <- "species"
zetaAnalysis <- read.table(paste("coastalBMIs",taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Calculate the modeled CS.
zetaModel <- lm(ConditionScore~zeta_1+zeta_2+zeta_10,data=zetaAnalysis)
zetaAnalysis$modeledConditionScore <- zetaModel$fitted.values
#Find common environmental parameters, with a VIF < 10, which contribute to the variation in the modeled CS.
environmentNames <- c(chemNames,sedimentNames,"meanDist","meanDepth","meanOutfallDist")
zetaModel <- lm(formula(paste("modeledConditionScore ~ ",paste(environmentNames,collapse="+"))),data=zetaAnalysis)
varList <- car::vif(zetaModel)
varList <- names(varList[varList <= 10])
zetaCor <- zetaAnalysis[,c("ConditionScore","modeledConditionScore","zeta_1","zeta_2","zeta_10",varList)]
zetaCor <- rcorr(as.matrix(zetaCor),type="pearson")
corr <- zetaCor$r
p.mat <- zetaCor$P
colnames(corr) <- c("CS","Modeled CS",":zeta[1]",":zeta[2]",":zeta[10]",varList)
colnames(corr)[colnames(corr)=="meanDepth"] <- "Depth"
colnames(corr)[colnames(corr)=="meanDist"] <- "Distance"
colnames(corr)[colnames(corr)=="meanOutfallDist"] <- "Outfall"
rownames(corr) <- c("CS","Modeled CS",":zeta[1]",":zeta[2]",":zeta[10]",varList)
rownames(corr)[rownames(corr)=="meanDepth"] <- "Depth"
rownames(corr)[rownames(corr)=="meanDist"] <- "Distance"
rownames(corr)[rownames(corr)=="meanOutfallDist"] <- "Outfall"
corrplot(corr = corr, p.mat = p.mat, diag = FALSE, type="upper", sig.level = 0.0001, tl.col="black", tl.srt=45, tl.cex=1.3, order="original", title=paste("BMIs aggregated to\n",taxonomicLevel),mar=c(0,0,3,0))

#Color coded scatterplots of zeta diversity analysis variables.
require(ggplot2)
require(viridis)
zetaPlot <- ggplot(zetaAnalysis, aes(x=ConditionScore,y=modeledConditionScore,color=PLAIC-ExpAIC))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=modeledConditionScore))
zetaPlot+xlab("Condition Score")+ylab("Modeled Condition Score")+scale_color_gradientn("Mean Depth (m)",colours = rev(plasma(10)))

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
#Prepare data input for mapping.
mapInput <- communityInput[,colnames(communityInput) %in% c("latitude","longitude","Stratum","DJG_Stratum1","Region","StationWaterDepth","ConditionScore",chemNames,sedimentNames)]
mapInput <- mapInput[!duplicated(mapInput),]
#Map data for continuous variables.
CalMap = leaflet(mapInput) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=viridis(10),domain=mapInput$Total_DDT)
CalMap %>% addCircleMarkers(color = ~ColorScale(Total_DDT), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~Total_DDT,title="SCCWRP Sample<br>Total DDT (ppm)")

#Map data for categorical variables.
pal <- colorFactor(palette = 'viridis', domain = mapInput$Stratum)
CalMap = leaflet(stationInput) %>% 
  addTiles()
CalMap %>% addCircleMarkers(color = ~pal(Stratum), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=pal,values=~Stratum,title="SCCWRP sample<br>Habitats")
