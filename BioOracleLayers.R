##This script extracts long term average environmental map layers from the BioOracle
##data sets and then extracts their values at particular sampling locations.
rm(list=ls())
require(dplyr)
require(tidyr)
require(zetadiv)
require(geosphere)
require(relaimpo)
require(tidyverse)
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
require(sdmpredictors)

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

#Filter out sample locations before 2000 or after 2014.
communityInput <- communityInput[communityInput$SampleYear >= 2000 & communityInput$SampleYear <= 2014,]

#Get all marine BioOracle layers.
env_data <- list_datasets(terrestrial = F,marine=T)
env_layers <- list_layers(env_data)

#Only keep data collected between 2000 and 2014.
env_layers <- env_layers[env_layers$end_year >= 2014 & env_layers$start_year <= 2000,]

#Only keep data containing long term average values at the sea surface or mean seafloor depth.
meanTerms <- c("mean value at mean bottom depth","mean value at sea surface")
env_layers <- subset(env_layers, derivation %in% meanTerms)

#Get SCCWRP sampling locations for 2000 - 2014.
mapInput <- communityInput[,c("StationID","latitude","longitude")]
mapInput <- mapInput[!duplicated(mapInput),]

#Reprojet sampling points into the same map projection as the environmental layers.
map_points <- SpatialPoints(mapInput[,c("longitude","latitude")],CRS("+proj=longlat +datum=WGS84"))
map_points <- spTransform(map_points, equalareaproj)

#Get map extent of SCCWRP sampling points.  Use this to clip the BioOracle map layers.
map_extent <- bbox(map_points)

#Download all map layers to be used into the working directory.
env_rasters <- load_layers(env_layers,equalarea=T,rasterstack=F,datadir=wd)

#Loop through all of the map layers and extract their values at the SCCWRP sampling points.
#Merge all of these values with the sample IDs and their locations.
env_values <- mapInput
for(layer_name in env_rasters){
  tmp <- raster::crop(layer_name,map_extent)
  env_value <- as.data.frame(raster::extract(tmp,map_points))
  colnames(env_value) <- layer_name@data@names
  env_values <- cbind(env_values,env_value)
}

#Output BioOracle layer values at SCCWRP sampling points to a file for use elsewhere.
write.table(env_values,"BioOracleLayers.txt",quote=FALSE,sep="\t",row.names = FALSE)
