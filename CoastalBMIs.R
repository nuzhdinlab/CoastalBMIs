rm(list=ls())
require(dplyr)
require(tidyr)
require(zetadiv)
require(geosphere)
require(relaimpo)
require(tidyverse)
require(pracma)
require(data.table)
require(Rfast)

wd <- "/scratch/alsimons/CoastalBMIs" #On the cluster
wd <- "~/Desktop/CoastalBMIs/" #Local drive
setwd(wd)

#Read in taxa abundance by sample.
communityInput <- read.table("Master - Benthos ed12.csv", header=T, sep=",",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")

#Create unique ID per sample merging station ID with sample replicate.
communityInput <- unite(communityInput,SampleID,c("StationID","Replicate"),sep="_",remove=FALSE)

#Create basic sample metadata names list.
sampleNames <- c("SampleID","StationWaterDepth","Outfall","latitude","longitude")

#Remove empty taxon rows.
communityInput <- communityInput[communityInput$Taxon!="",]

#Remove duplicate entries.
communityInput <- communityInput[!duplicated(communityInput[,c("SampleID","Taxon")]),]

#Read in sample station metadata.
stationInput <- read.table("Master Data - Station Info.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8",allowEscapes=TRUE)

#Merge in sample sation metadata to community data.
communityInput <- dplyr::left_join(communityInput,stationInput,by=c("StationID"))

#Read in distance to outfall pipe per sample location.
outfallInput <- read.table("SampleStationsWithDistanceToOutfall.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8",allowEscapes=TRUE)
#Rename ID field
outfallInput <- dplyr::rename(outfallInput, StationID := ID)
#Calculate the sum of the error function of the inverse square root of the distances between each outfall pipe
#and a given sampling location.  Use this to estimate an outfall concentration at each sampling point.
varList <- colnames(outfallInput[,2:ncol(outfallInput)])
outfallInput[,varList] <- lapply(outfallInput[,varList], function(x) erf(1/sqrt(x)))
outfallInput$Outfall <- rowSums(outfallInput[,varList])
#Keep only relevant columns.
outfallInput <- outfallInput[,c("StationID","Outfall")]

#Merge in distance to outfall pipe locations into the community data.
communityInput <- dplyr::left_join(communityInput,outfallInput,by=c("StationID"))

#Read in BMI taxonomic data.
taxaInput <- read.table("Master - Taxonomic Info edit.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8",quote="")
#Remove quote characters.
taxaInput <- as.data.frame(sapply(taxaInput, function(x) gsub("\"", "", x)))
#Coerce Taxon column to character
taxaInput$Taxon <- as.character(taxaInput$Taxon)

#Merge in taxonomic data into community data.
communityInput <- dplyr::left_join(communityInput,taxaInput,by=c("Taxon"))

#Get the taxonomic richness per sample at various levels of taxonomic aggregation.
setDT(communityInput)[, Nspecies:=uniqueN(species), by=SampleID]
setDT(communityInput)[, NGenera:=uniqueN(Genus), by=SampleID]
setDT(communityInput)[, NFamilies:=uniqueN(Family), by=SampleID]
setDT(communityInput)[, NOrders:=uniqueN(Order), by=SampleID]
setDT(communityInput)[, NClasses:=uniqueN(Class), by=SampleID]
setDT(communityInput)[, NPhyla:=uniqueN(Phylum), by=SampleID]

#Read in assessment scores.
assessmentInput <- read.table("Master Data - Assessment Scores.csv", header=T, sep=",", as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Create unique ID per sample merging station ID with sample replicate.
assessmentInput <- unite(assessmentInput,SampleID,c("StationID","Replicate"),sep="_",remove=FALSE)

#Keep only the benthic response index (BRI) or sediment quality objectives benthic response index (SQO_BRI)
assessmentInput <- assessmentInput[assessmentInput$Index %in% c("BRI","SQO_BRI"),]

#Reformat assessment data so the each unique assessment score is a column.
assessmentInput <- assessmentInput %>% dplyr::group_by(Index) %>% dplyr::mutate(n = 1:n()) %>% tidyr::spread(Index,Index_Score)
assessmentInput$n <- NULL
#Get assessment names.
assessmentNames <- c("BRI","SQO_BRI")

#Keep only unique Condition Scores per sample.
assessmentInput <- assessmentInput[,c("SampleID","ConditionScore",assessmentNames)]
assessmentInput <- assessmentInput[!duplicated(assessmentInput),]

#Merge in assessment scores into community data.
communityInput <- dplyr::left_join(communityInput,assessmentInput,by=c("SampleID"))

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
taxonomicLevel <- c("Phylum") #Choose a taxonomic level to aggregate count data on.

#Create simplified community data frame to summarize data at a particular taxonomic level.
communityInputSummarized <- dplyr::select(communityInput,SampleID,taxonomicLevel,Abundance)
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

#tmp <- communityInput[,c("SampleID","SQO_BRI","BRI",chemNames,"Outfall","Nspecies","NGenera","NFamilies","NOrders","NClasses","NPhyla")]
#tmp <- tmp[!duplicated(tmp),]
#write.table(tmp,"OutfallAndCommunityData.csv",quote=FALSE,sep=",",row.names = FALSE)

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
          metadataSubset <- metadataSubset[,colnames(metadataSubset) %in% c(sampleNames,chemNames,sedimentNames,assessmentNames)]
          metadataSubset <- metadataSubset[!duplicated(metadataSubset),]
          meanDist <- mean(distm(metadataSubset[,c("longitude","latitude")]))
          sdDist <- sd(distm(metadataSubset[,c("longitude","latitude")]))
          zeta_2_DD <- Zeta.ddecay(xy=metadataSubset[,c("longitude","latitude")],data.spec=data.spec,order=2,distance.type="ortho",rescale=T,normalize="Jaccard",trsf="log",plot=F)
          zeta_2_DD <- zeta_2_DD$reg$coefficients[2]
          zeta_10_DD <- Zeta.ddecay(xy=metadataSubset[,c("longitude","latitude")],data.spec=data.spec,order=10,distance.type="ortho",rescale=T,normalize="Jaccard",trsf="log",plot=F)
          zeta_10_DD <- zeta_10_DD$reg$coefficients[2]
          meanOutfall <- mean(metadataSubset$Outfall)
          sdOutfall <- sd(metadataSubset$Outfall)
          meanDepth <- mean(metadataSubset$StationWaterDepth)
          sdDepth <- sd(metadataSubset$StationWaterDepth)
          chemData <- metadataSubset[,colnames(metadataSubset) %in% chemNames]
          chemRow <- colMeans(chemData,na.rm=T)
          chemRow[is.nan(chemRow)] <- NA
          sedimentData <- metadataSubset[,colnames(metadataSubset) %in% sedimentNames]
          sedimentRow <- colMeans(sedimentData,na.rm=T)
          sedimentRow[is.nan(sedimentRow)] <- NA
          assessmentData <- metadataSubset[,colnames(metadataSubset) %in% assessmentNames]
          assessmentRow <- colMeans(assessmentData,na.rm=T)
          assessmentRow[is.nan(assessmentRow)] <- NA
          envVar <- colVars(as.matrix(cbind(chemData,sedimentData,metadataSubset[,c("StationWaterDepth","Outfall")])),na.rm=T)
          envVar <- sqrt(sum(envVar[!is.nan(envVar)]^2))/mean(envVar[!is.nan(envVar)])
          dataRow <- t(cbind(t(as.data.frame(list(c(stratumType,year,score,meanDepth,sdDepth,meanDist,sdDist,zeta_2_DD,zeta_10_DD,meanOutfall,sdOutfall,zetaOrders,zetaSDOrders,ExpExp,ExpAIC,PLExp,PLAIC,as.list(as.numeric(sedimentRow)),as.list(as.numeric(chemRow)),as.list(as.numeric(assessmentRow)),envVar))))))
          colnames(dataRow) <- c("stratumType","year","ConditionScore","meanDepth","sdDepth","meanDist","sdDist","zeta_2_DD","zeta_10_DD","meanOutfall","sdOutfall",zetaNames,zetaSDNames,"ExpExp","ExpAIC","PLExp","PLAIC",sedimentNames,chemNames,assessmentNames,"envVar")
          rownames(dataRow) <- NULL
          zetaAnalysis <- rbind(zetaAnalysis,dataRow)
          print(paste(j,stratumType,year,score,meanDepth,sdDepth,meanDist,sdDist,zeta_2_DD,zeta_10_DD,meanOutfall,sdOutfall,zetaOrders,zetaSDOrders,ExpExp,ExpAIC,PLExp,PLAIC,envVar))
          print(chemRow)
          print(sedimentRow)
          print(assessmentRow)
        }
      }
    }
  }
}

colnames(zetaAnalysis) <- c("stratumType","year","ConditionScore","meanDepth","sdDepth","meanDist","sdDist","zeta_2_DD","zeta_10_DD","meanOutfall","sdOutfall",zetaNames,zetaSDNames,"ExpExp","ExpAIC","PLExp","PLAIC",sedimentNames,chemNames,assessmentNames,"envVar")
tmp <- zetaAnalysis$stratumType
indx <- sapply(zetaAnalysis, is.factor)
zetaAnalysis[indx] <- lapply(zetaAnalysis[indx], function(x) as.numeric(as.character(x)))
zetaAnalysis <- do.call(data.frame,lapply(zetaAnalysis, function(x) replace(x, is.infinite(x),NA)))
zetaAnalysis$stratumType <- tmp

#Calculate the modeled SQO BRI scores.
tmp1 <- zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),]
zetaModel <- lm(SQO_BRI~zeta_1+zeta_2+zeta_10,data=tmp1)
tmp1$modeledSQO_BRI <- zetaModel$fitted.values
tmp1$AssessmentType <- "SQO_BRI"
tmp1$modeledBRI <- NA
#Calculate the modeled BRI scores.
tmp2 <- zetaAnalysis[!is.na(zetaAnalysis$BRI),]
zetaModel <- lm(BRI~zeta_1+zeta_2+zeta_10,data=tmp2)
tmp2$modeledBRI <- zetaModel$fitted.values
tmp2$AssessmentType <- "BRI"
tmp2$modeledSQO_BRI <- NA

#Merge in data frames together for zeta diversity analysis.
zetaAnalysis <- rbind(tmp1,tmp2)
#Save zeta diversity analysis for a given taxonomic level.
write.table(zetaAnalysis,paste("coastalBMIs",taxonomicLevel,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#

##To run locally.
taxonomicLevel <- "Genus"
zetaAnalysis <- read.table(paste("coastalBMIs",taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Determine the relative importance of various orders of zeta diversity in models of biotic integrity.
calc.relimp(lm(SQO_BRI~zeta_1+zeta_2+zeta_10,data=zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),]))
calc.relimp(lm(BRI~zeta_1+zeta_2+zeta_10,data=zetaAnalysis[!is.na(zetaAnalysis$BRI),]))

#Assess contributions of environmental variables to modeled measures of biotic integrity.
require(car)
environmentNames <- c("meanDist","meanDepth","meanOutfall","envVar",sedimentNames,chemNames)
zetaModel <- lm(formula(paste("modeledSQO_BRI ~ ",paste(environmentNames,collapse="+"))),data=zetaAnalysis)
calc.relimp(zetaModel)
#Filter our variables with a high variance inflation factor.
varList <- car::vif(zetaModel)
varList <- names(varList[varList <= 10])
#Model the modeled condition score as a function of the remaining environmental metrics
zetaModel <- lm(formula(paste("modeledSQO_BRI ~ ",paste(varList,collapse="+"))),data=zetaAnalysis)
calc.relimp(zetaModel)

#Assess contributions of environmental variables to modeled measures of biotic integrity.
require(car)
environmentNames <- c(chemNames,sedimentNames,"meanDist","meanDepth","meanOutfall")
zetaModel <- lm(formula(paste("SQO_BRI ~ ",paste(environmentNames,collapse="+"))),data=zetaAnalysis)
calc.relimp(zetaModel)
#Filter our variables with a high variance inflation factor.
varList <- car::vif(zetaModel)
varList <- names(varList[varList <= 10])
#Model the modeled condition score as a function of the remaining environmental metrics
zetaModel <- lm(formula(paste("SQO_BRI ~ ",paste(varList,collapse="+"))),data=tmp)
calc.relimp(zetaModel)

# Assessing R2 shrinkage using 10-Fold Cross-Validation 
require(bootstrap)
require(caret)
set.seed(1)
train.control <- trainControl(method="repeatedcv",number=10,repeats=10)
ConditionScoremodel <- train(SQO_BRI~zeta_1+zeta_2+zeta_10,data=zetaAnalysis,na.action=na.omit,method="lm",trControl=train.control)
print(ConditionScoremodel)
ConditionScoremodel <- train(BRI~zeta_1+zeta_2+zeta_10,data=zetaAnalysis,na.action=na.omit,method="lm",trControl=train.control)
print(ConditionScoremodel)

#Individual correlation plots of mean and modeled condition scores against environmental parameters
#and various measures of zeta diversity.
require(Hmisc)
require(corrplot)
#Find common environmental parameters which contribute to the variation in the modeled CS.
environmentNames <- c("meanDepth","meanOutfall","envVar")
#Create correlation matrix, with significance values, for SQO_BRI scores.
zetaCor <- zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),c("SQO_BRI","modeledSQO_BRI","zeta_1","zeta_2","zeta_10","PLAIC","ExpAIC",environmentNames)]
#Calculate PLAIC - ExpAIC, a measure of the relative likelihood of niche differentiation versus stochastic community assembly processes.
zetaCor$deltaAIC <- (zetaCor$PLAIC - zetaCor$ExpAIC)
zetaCor$PLAIC <- NULL
zetaCor$ExpAIC <- NULL
zetaCor <- rcorr(as.matrix(zetaCor),type="spearman")
corr <- zetaCor$r
p.mat <- zetaCor$P
#Rename variables for plotting.
corr <- as.data.frame(corr)
setnames(corr, old=c("meanDepth","meanOutfall","SQO_BRI","modeledSQO_BRI","zeta_1","zeta_2","zeta_10","deltaAIC"), new=c("Depth","Outfall","SQO BRI","Modeled SQO BRI",":zeta[1]",":zeta[2]",":zeta[10]",":Delta[AIC]"))
corr <- as.matrix(corr)
#Create correlation plot.
corrplot(corr = corr, p.mat = p.mat, diag = FALSE, type="upper",
         sig.level = 0.0001, tl.col="black", tl.srt=45, tl.cex=1.3,
         order="original", title=paste("BMIs with SQO BRI scores\naggregated to",taxonomicLevel),
         mar=c(0,0,3,0))

#Create correlation matrix, with significance values, using BRI scores.
zetaCor <- zetaAnalysis[!is.na(zetaAnalysis$BRI),c("BRI","modeledBRI","zeta_1","zeta_2","zeta_10","PLAIC","ExpAIC",environmentNames)]
#Calculate PLAIC - ExpAIC, a measure of the relative likelihood of niche differentiation versus stochastic community assembly processes.
zetaCor$deltaAIC <- (zetaCor$PLAIC - zetaCor$ExpAIC)
zetaCor$PLAIC <- NULL
zetaCor$ExpAIC <- NULL
zetaCor <- rcorr(as.matrix(zetaCor),type="spearman")
corr <- zetaCor$r
p.mat <- zetaCor$P
#Rename variables for plotting.
corr <- as.data.frame(corr)
setnames(corr, old=c("meanDepth","meanOutfall","BRI","modeledBRI","zeta_1","zeta_2","zeta_10","deltaAIC"), new=c("Depth","Outfall","BRI","Modeled BRI",":zeta[1]",":zeta[2]",":zeta[10]",":Delta[AIC]"))
corr <- as.matrix(corr)
#Create correlation plot.
corrplot(corr = corr, p.mat = p.mat, diag = FALSE, type="upper",
         sig.level = 0.0001, tl.col="black", tl.srt=45, tl.cex=1.3,
         order="original", title=paste("BMIs with BRI scores\naggregated to",taxonomicLevel),
         mar=c(0,0,3,0))

#Multipanel correlation plots of mean and modeled condition scores against environmental parameters
#and various measures of zeta diversity.
#Create a merged zeta diversity data frame for all taxonomic levels of aggregation.
zetaTotal <- data.frame()
for(taxonomicLevel in taxonomicLevels){
  zetaAnalysis <- read.table(paste("coastalBMIs",taxonomicLevel,".txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
  zetaAnalysis$Level <- taxonomicLevel
  zetaAnalysis$deltaAIC <- (zetaAnalysis$PLAIC - zetaAnalysis$ExpAIC)
  zetaTotal <- rbind(zetaAnalysis,zetaTotal)
}
zetaTotal$PLAIC <- NULL
zetaTotal$ExpAIC <- NULL

#Split nearshore zeta diversity data into a list of data frames using taxonomic level as the split variable. 
B <- split(zetaTotal[!is.na(zetaTotal$SQO_BRI),c("SQO_BRI","modeledSQO_BRI","zeta_1","zeta_2","zeta_10","deltaAIC",environmentNames)],zetaTotal[!is.na(zetaTotal$SQO_BRI),"Level"])
M <- lapply(B, function(x) rcorr(as.matrix(x),type="spearman"))
M <- M[taxonomicLevels]

#Create multipanel nearshore zeta diversity correlation plots.
dev.off()
par(mfrow=c(3,2))
col<- colorRampPalette(c("red","white","blue"))(40)
i=1
for(Level in names(M)){
  corr <- M[[Level]]$r
  p.mat <- M[[Level]]$P
  #Rename variables for plotting.
  corr <- as.data.frame(corr)
  setnames(corr, old=c("meanDepth","meanOutfall","SQO_BRI","modeledSQO_BRI","zeta_1","zeta_2","zeta_10","deltaAIC"), new=c("Depth","Outfall","SQO BRI","Modeled SQO BRI",":zeta[1]",":zeta[2]",":zeta[10]",":Delta[AIC]"))
  corr <- as.matrix(corr)
  #Create correlation plot.
  corrplot(corr = corr, p.mat = p.mat, diag = FALSE, type="lower",
           sig.level = 0.0001, tl.col="black", tl.srt=15, tl.cex=0.8,
           order="original",mar=c(0,0,1,0))
  mtext(paste("(",LETTERS[i],")",sep=""), side=3, line=1, cex=2, adj=0, las=0)
  mtext(Level, side=2, line=0, cex=1.5, adj=0.5, las=3)
  box("figure", col="black", lwd = 2)
  i=i+1
}

#Split offhore zeta diversity data into a list of data frames using taxonomic level as the split variable. 
B <- split(zetaTotal[!is.na(zetaTotal$BRI),c("BRI","modeledBRI","zeta_1","zeta_2","zeta_10","deltaAIC",environmentNames)],zetaTotal[!is.na(zetaTotal$BRI),"Level"])
M <- lapply(B, function(x) rcorr(as.matrix(x),type="spearman"))
M <- M[taxonomicLevels]

#Create multipanel offshore zeta diversity correlation plots.
dev.off()
par(mfrow=c(3,2))
col<- colorRampPalette(c("red","white","blue"))(40)
i=1
for(Level in names(M)){
  corr <- M[[Level]]$r
  p.mat <- M[[Level]]$P
  #Rename variables for plotting.
  corr <- as.data.frame(corr)
  setnames(corr, old=c("meanDepth","meanOutfall","BRI","modeledBRI","zeta_1","zeta_2","zeta_10","deltaAIC"), new=c("Depth","Outfall","BRI","Modeled BRI",":zeta[1]",":zeta[2]",":zeta[10]",":Delta[AIC]"))
  corr <- as.matrix(corr)
  #Create correlation plot.
  corrplot(corr = corr, p.mat = p.mat, diag = FALSE, type="lower",
           sig.level = 0.0001, tl.col="black", tl.srt=15, tl.cex=0.8,
           order="original",mar=c(0,0,1,0))
  mtext(paste("(",LETTERS[i],")",sep=""), side=3, line=1, cex=2, adj=0, las=0)
  mtext(Level, side=2, line=0, cex=1.5, adj=0.5, las=3)
  box("figure", col="black", lwd = 2)
  i=i+1
}

#Check contributions to variations in how zeta diversity decays with distance.
summary(aov(zeta_2_DD~stratumType*year*SQO_BRI,data=zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),]))
summary(aov(zeta_10_DD~stratumType*year*SQO_BRI,data=zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),]))
summary(aov(zeta_2_DD~stratumType*year*BRI,data=zetaAnalysis[!is.na(zetaAnalysis$BRI),]))
summary(aov(zeta_10_DD~stratumType*year*BRI,data=zetaAnalysis[!is.na(zetaAnalysis$BRI),]))

#Compare steepness of zeta diversity decay with distance between nearshore and offshore environments.
wilcox.test(zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),"zeta_2_DD"],zetaAnalysis[!is.na(zetaAnalysis$BRI),"zeta_2_DD"],alternative="g")
wilcox.test(zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),"zeta_10_DD"],zetaAnalysis[!is.na(zetaAnalysis$BRI),"zeta_10_DD"],alternative="g")

#Compare environmental variability between nearshore and offshore environments.
wilcox.test(zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),"envVar"],zetaAnalysis[!is.na(zetaAnalysis$BRI),"envVar"],alternative="t")

#Relative likelihoods of community assembly models
wilcox.test(zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),"PLAIC"],zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),"ExpAIC"],alternative="l")
wilcox.test(zetaAnalysis[!is.na(zetaAnalysis$BRI),"PLAIC"],zetaAnalysis[!is.na(zetaAnalysis$BRI),"ExpAIC"],alternative="l")
t.test(zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),"PLAIC"],zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),"ExpAIC"],alternative="l")
t.test(zetaAnalysis[!is.na(zetaAnalysis$BRI),"PLAIC"],zetaAnalysis[!is.na(zetaAnalysis$BRI),"ExpAIC"],alternative="l")

#Comparing levels of environmental heterogeneity.
wilcox.test(zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),"envVar"],zetaAnalysis[!is.na(zetaAnalysis$BRI),"envVar"])

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
mapInput <- communityInput[,colnames(communityInput) %in% c("SampleID","latitude","longitude","SampleYear","Stratum","DJG_Stratum1","Region","StationWaterDepth","ConditionScore","Outfall",chemNames,sedimentNames,assessmentNames)]
mapInput <- mapInput[mapInput$SampleID %in% sampleList$SampleID,]
mapInput <- mapInput[!duplicated(mapInput),]

mapInput <- mapInput[!is.na(mapInput$SQO_BRI),]
#Map data for continuous variables.
CalMap = leaflet(mapInput) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=plasma(10),domain=mapInput$SQO_BRI)
CalMap %>% addCircleMarkers(color = ~ColorScale(SQO_BRI), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~SQO_BRI,title="SQO BRI")

#Map data for categorical variables.
pal <- colorFactor(palette = 'viridis', domain = mapInput$Stratum)
CalMap = leaflet(mapInput) %>% 
  addTiles()
CalMap %>% addCircleMarkers(color = ~pal(Stratum), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=pal,values=~Stratum,title="SCCWRP sample<br>habitat types")
