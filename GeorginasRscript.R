getwd()
install.packages("stats")
require(stats)
my_data <- read.delim("coastalBMIsPhylum.txt")
ZetaAnalysis <- read.delim("coastalBMIsPhylum.txt")
my_data
plot(ZetaAnalysis)
ZetaAnalysis
plot(ConditionScore, zeta_1)
data('coastalBMIsPhylum')
data('ConditionScore')


plot(x=ConditionScore, y=zeta_1) 

zetaAnalysis <- read.delim('coastalBMIsPhylum.txt', sep = '\t', header = T)
plot(zetaAnalysis$ConditionScore, zetaAnalysis$zeta_1)

plot(zetaAnalysis$ConditionScore, zetaAnalysis$zeta_1,
     xlab="Condition Score",
     ylab="zeta_1")

aov(y ~ A + B + A:B, data=mydataframe)

a = aov(ConditionScore~zeta_1, data=zetaAnalysis)
summary(a)


a = aov(zeta_1~ConditionScore, data=zetaAnalysis)
summary(a)

zetaAnalysis <- read.delim("coastalBMIsPhylum.txt")
!is.na(zetaAnalysis$SQO_BRI)
zetaAnalysis<- zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),]
nrow(zetaAnalysis)


a = aov(zeta_2_DD~stratumType*year*SQO_BRI+stratumType:year:SQO_BRI, data=zetaAnalysis)
summary(a)

citation("zetadiv")
install.packages("zetadiv")


summary(aov(zeta_2_DD~stratumType*year*SQO_BRI,data=zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),]))
summary(aov(zeta_10_DD~stratumType*year*SQO_BRI,data=zetaAnalysis[!is.na(zetaAnalysis$SQO_BRI),])) 

summary(aov(zeta_2_DD~stratumType*year*BRI,data=zetaAnalysis[!is.na(zetaAnalysis$BRI),]))
summary(aov(zeta_10_DD~stratumType*year*BRI,data=zetaAnalysis[!is.na(zetaAnalysis$BRI),])) 
zetaAnalysis <- read.delim("coastalBMIsClass.txt")
zetaAnalysis <- read.delim("coastalBMIsOrder.txt")
zetaAnalysis <- read.delim("coastalBMIsFamily.txt")
zetaAnalysis <- read.delim("coastalBMIsGenus.txt")
zetaAnalysis <- read.delim("coastalBMIsspecies.txt")



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


setwd("/Users/georginagemayel/Desktop/SummerR")

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

#Get the percent occurrence of each sediment type.
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
8mapInput <- communityInput[,colnames(communityInput) %in% c("SampleID","latitude","longitude","SampleYear","Stratum","DJG_Stratum1","Region","StationWaterDepth","ConditionScore","Outfall",chemNames,sedimentNames,assessmentNames)]
mapInput <- mapInput[mapInput$SampleID %in% mapInput$SampleID,]
mapInput <- mapInput[!duplicated(mapInput),]
mapInput <- mapInput[!is.na(mapInput$SQO_BRI),]
#Map data for continuous variables.
CalMap = leaflet(mapInput) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=plasma(10),domain=mapInput$SQO_BRI)
CalMap %>% addCircleMarkers(color = ~ColorScale(SQO_BRI), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~SQO_BRI,title="SQO BRI")

#try
mapInput <- communityInput[,colnames(communityInput) %in% c("SampleID","latitude","longitude","SampleYear","Stratum","DJG_Stratum1","Region","StationWaterDepth","ConditionScore","Outfall",chemNames,sedimentNames,assessmentNames)]
mapInput <- mapInput[mapInput$SampleID %in% mapInput$SampleID,]
mapInput <- mapInput[!duplicated(mapInput),]
mapInput <- mapInput[!is.na(mapInput$BRI),]

CalMap = leaflet(mapInput) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=plasma(10),domain=mapInput$BRI)
CalMap %>% addCircleMarkers(color = ~ColorScale(BRI), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~BRI,title="BRI")

#Map data for categorical variables.
pal <- colorFactor(palette = 'viridis', domain = mapInput$Stratum)
CalMap = leaflet(mapInput) %>% 
  addTiles()
CalMap %>% addCircleMarkers(color = ~pal(Stratum), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=pal,values=~Stratum,title="SCCWRP sample<br>habitat types")

#try

pal <- colorFactor(palette = 'viridis', domain = mapInput$SampleYear)
CalMap = leaflet(mapInput) %>% 
  addTiles()
CalMap %>% addCircleMarkers(color = ~pal(SampleYear), fill = TRUE,radius=1,fillOpacity = 1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=pal,values=~SampleYear,title="SCCWRP sample<br>year")


taxonomicLevels <- c("Phylum","Class","Order","Family","Genus","species")
taxonomicLevel <- c("Phylum") #Choose a taxonomic level to aggregate count data on.


nrow(communityInput)

colnames(communityInput)
Phylum_column <- communityInput$Phylum
Phylum_column
length(Phylum_column)
Phylum_unique <- unique(Phylum_column)
Phylum_unique
length(Phylum_unique)

#use 
length(unique(communityInput$Phylum))

left_join(mapInput, ZetaAnalysis, by= c("SampleYear" = "year"))
left_join(mapInput, ZetaAnalysis, by= c("Stratum" = "stratumType"))

mapInputDistance <- left_join(mapInput, ZetaAnalysis, by= c("SampleYear" = "year", "Stratum" = "stratumType"))
